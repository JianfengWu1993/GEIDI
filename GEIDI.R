library(methods)

options(stringsAsFactors = F)

args <- commandArgs()

infile = "D:\\ADNI_chow\\hippotest.txt"

gfile = "D:\\ADNI_chow\\group.txt"

cfile = "D:\\ADNI_chow\\casetest.txt" # casefile #

outfile = "D:\\ADNI_chow\\chow_result.txt"

#outdir = args[10] #"D:\\ADNI_Gene\\Chow_test\\"# output dir #

pcut = 1

dofdr = TRUE

#plotchow = {{args.plot | R}}

#devpars  = {{args.devpars | R}}

#ggs      = {{args.ggs | R}}

#inopts   = {{args.inopts | R}}

covfile = ''

if (dofdr == T) dofdr = 'BH'

bQuote = function(s) {
    paste0('`', s, '`')
}

regress = function(regdata, name, fmula = NULL) {

	Y = bQuote(colnames(regdata)[ncol(regdata)])

	fmula = ifelse(is.null(fmula), paste(Y, "~ ."), fmula)
	
	fmula = unlist(strsplit(fmula, "~"))
	
	fmula = paste0(rev(fmula), collapse = "~")

	m = lm(as.formula(fmula), data = regdata)
  
	list(model = m, ssr = sum(m$residuals ^ 2), n = nrow(regdata), name = name, p_val = coef(summary(m))[2,4])

}



formlm = function(model, k, withname = T) {

	lencoef = length(model$model$coefficients)

	coefns = c(names(model$model$coefficients)[(lencoef - k + 2):lencoef], '_')

	coeffs = as.numeric(c(model$model$coefficients[(lencoef - k + 2):lencoef], model$model$coefficients[1]))

	coefns = c(coefns, "N")

	coeffs = c(coeffs, model$n)

	if (withname) {

		paste0(model$name, ': ', paste(coefns, round(coeffs, 3), sep = '=', collapse = ' '))

	} else {

		paste(coefns, round(coeffs, 3), sep = '=', collapse = ' ')

	}

}



chow = function(pooled, subregs, k, case) {

	subssr = sum(sapply(subregs, function(x) x$ssr))

	ngroups = length(subregs)

	J = (ngroups - 1) * k

	DF = pooled$n - ngroups * k

	FS = (pooled$ssr - subssr) * DF / subssr / J

	groups = lapply(subregs, function(m) formlm(m, k))

	pooledm = formlm(pooled, k, FALSE)

	list(Case = case, Pooled = pooledm, Groups = paste(groups, collapse = '; '), all_pval = pooled$p_val, group_pval = paste(sapply(subregs, function(x) x$p_val), collapse = '; '), Fstat = FS, Pval = pf(FS, J, DF, lower.tail = FALSE))

}



model2eq = function(model) {

	vars = colnames(model$model)

	cf = sapply(model$coefficients, function(f) {

		if (is.na(f)) return("NA")

		f = round(f, 2)

		if (f >= 0) return(paste0("+", f))

		return(paste0("-", - f))

	})

	paste0(vars[1], ' = ', cf[1], paste(

		sapply(

			2:length(cf),

			function(x) paste0(cf[x], '*', vars[x])

		), collapse = ''

	))

}

# format data.frame to output
pretty.numbers = function(df, formats) {
    # remember set stringsAsFactors as FALSE for the dataframe!!
    if (nrow(df) == 0) {
        return(df)
    }
    allCols = colnames(df)
    formatedCols = c()
    for (fcols in names(formats)) {
        if (fcols == '.') {
            # must be last element of formats
            cols = which(!allCols %in% formatedCols)
        } else {
            cols = unlist(strsplit(fcols, '..', fixed = T))
            formatedCols = c(formatedCols, cols)
        }
        cols = intersect(cols, allCols)
        df[, cols] = sprintf(formats[[fcols]], as.numeric(unlist(df[, cols])))
    }
    df
}

#indata = read.table.inopts(infile, inopts)

indata = read.table(infile, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

#     X1  X2  X3  X4 ... Y

# G1  1   2   1   4  ... 9

# G2  2   3   1   1  ... 3

# ... ...

# Gm  3   9   1   7  ... 8

#K = ncol(indata)

if (covfile != "") {

	covdata = read.table(covfile, header = T, row.names = 1, check.names = F)

	indata = cbind(covdata[rownames(indata),, drop = F], indata)

}

gdata = read.table(gfile, header = F, row.names = NULL, check.names = F)

# G1    Group1  Case1

# G2    Group1  Case1

# ... ...

# Gs    Group2  Case1

# Gt    Group2  Case1

# Gt    Group1  Case2

# ... ...

# Gu    Group1  Case2

# ... ...

# Gz    Group2  Case2
cases = ''
if (ncol(gdata) == 2) {
	# no case

	cases = 'Case1'

}

fmulas = NULL

if (cfile != "") {

	fmulas = read.table(cfile, header = F, row.names = 1, sep = "\t", check.names = F)

	cases = rownames(fmulas)

}



results = data.frame(

	Case = character(),

	Pooled = character(),

	Groups = character(),
	
	pooled_pval = double(),
	
	group_pval = character(),
	
	Fstat = double(),

	Pval = double()
	

	)


for (case in cases) {

	caserows = if (ncol(gdata) > 2) gdata[which(gdata$Case == case),, drop = F] else gdata

	fmula = fmulas[case,]
	
	allvars = all.vars(as.formula(fmula))

	K = ifelse(is.null(fmula), ncol(indata), length(allvars))

	pooled_lm = regress(indata[caserows[, 1],, drop = F], name = 'Pooled', fmula = fmulas[case,])

	subgroups = levels(factor(caserows[, 2]))

	subgrp_lm = lapply(subgroups, function(sgroup) {

		sdata = indata[caserows[which(caserows[, 2] == sgroup), 1],, drop = F]

		if (nrow(sdata) < 3) NULL

		else regress(sdata, name = sgroup, fmula = fmulas[case,])

	})

	subgrp_lm[sapply(subgrp_lm, is.null)] <- NULL

	# no subgroups
  
	if (length(subgrp_lm) < 2) next

	ret = chow(pooled_lm, subgrp_lm, k = K, case = case)

	if (is.na(ret$Pval) || ret$Pval >= pcut) next
  
	
	results = rbind(results, ret)
}

if (dofdr != F) {
    results = cbind(results, Qval = p.adjust(results$Pval, method = dofdr))
}
print(outfile)
print(results$Pval)
write.csv(pretty.numbers(results, list(
    Fstat = '%.3f',
    Pval..Qval = '%.3E'
)), file = outfile, row.names = F, quote = F)