rule xtailreport:
    input:
        fcplot="xtail/fc_{contrast}.pdf",
        rplot="xtail/r_{contrast}.pdf",
    output:
        fcplot=report("figures/fc_{contrast}.jpg", caption="../report/xtail_fc.rst", category="Regulation"),
        rplot=report("figures/r_{contrast}.jpg", caption="../report/xtail_r.rst", category="Regulation")
    conda:
        "../envs/imagemagick.yaml"
    threads: 1
    shell: ("mkdir -p figures; convert -density 150 -trim {input.fcplot}  -quality 100  -flatten -sharpen 0x1.0 {output.fcplot}; convert -density 150 -trim {input.rplot}  -quality 100  -flatten -sharpen 0x1.0 {output.rplot}; convert -density 150 -trim {input.fcplot}  -quality 100  -flatten -sharpen 0x1.0 {output.fcplot}; convert -density 150 -trim {input.rplot}  -quality 100  -flatten -sharpen 0x1.0 {output.rplot}; ")

rule ribotishreport:
    input:
        "ribotish/{condition}-{replicate}-qual.pdf" 
    output:
        report("figures/{condition}-{replicate}-qual.jpg", caption="../report/ribotishquality.rst", category="Ribotish")
    conda:
        "../envs/imagemagick.yaml"
    threads: 1
    shell: ("mkdir -p figures; convert -density 150 -trim {input} -quality 100 -flatten -sharpen 0x1.0 {output};")

