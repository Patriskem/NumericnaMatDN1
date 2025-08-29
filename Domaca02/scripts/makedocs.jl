using Weave

#Weave.weave("demo.jl", doctype="minted2pdf", out_path=:pwd)
Weave.weave("demo.jl", doctype="md2pdf", out_path=:pwd)