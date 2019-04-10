(TeX-add-style-hook
 "flagme"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("caption" "tableposition=top") ("inputenc" "utf8")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "amsmath"
    "amscd"
    "caption"
    "ifthen"
    "inputenc"))
 :latex)

