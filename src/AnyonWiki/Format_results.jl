#=----------------------------------------------------------
    Format the results of the computations around the
    AnyonWiki
----------------------------------------------------------=#

K = anyonwiki_keys(5)

s = IOBuffer()

write(s, """
\\begin{longtable}{c|c|c|p{5cm}}
\t\\textbf{Code} & \\textbf{Rank} & \\textbf{Multiplicity} & \\textbf{Splits over cyclotomic extension} \\\\ \\hline
""")

for k ∈ K
    meta = anyonwiki_center_meta(k...)
    F = meta["field"]
    split = F == QQ ? "Yes" : (is_cyclotomic_polynomial(F.pol) ? "Yes" : "No")
    write(s, "\t\\AnyonCode{$(join(k, ", "))} & $(meta["rank"]) & $(meta["multiplicity"]) & $(split) \\\\ \n")
end

s = String(take!(s)) * "\\end{longtable}"
