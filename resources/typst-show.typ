// Project override of the orange-book `typst-show.typ` partial (Quarto 1.9).
// Differences from the built-in partial:
//
// - Cover: the template's fixed-height title box overflows with a long title
//   and many authors. Its text is hidden (title/subtitle wrapped in `hide`,
//   the author string removed with a show rule) and its fill set to none, so
//   title and author are still passed on for the PDF metadata. The visible
//   cover is drawn through `cover:`, a full-page block with a coloured stripe
//   in the left margin and the title, authors and date to the right of it.
// - `copyright:` prints authors with numbered affiliations on the page after
//   the cover. Authors with `equal-contributor: true` are marked with a dagger
//   and a note, and the corresponding author with an asterisk.
#import "@preview/orange-book:0.7.1": book, part, chapter, appendices

#let cover-stripe-width = 25mm
#let cover-stripe-color = rgb("#b5cfb6")
// Accent colour for headings, the table of contents and links
#let accent-color = rgb("#056731")

#show: doc => {
$if(by-author)$
  show "$for(by-author)$$it.name.literal$$sep$, $endfor$": none
$endif$
  book.with(
$if(title)$
  title: [#hide[$title$]],
$endif$
$if(subtitle)$
  subtitle: [#hide[$subtitle$]],
$endif$
  cover-background: none,
  cover: block(width: 100%, height: 100%, inset: (left: cover-stripe-width + 15mm, right: 25mm))[
    #place(left + top, dx: -(cover-stripe-width + 15mm), rect(width: cover-stripe-width, height: 100%, fill: cover-stripe-color))
    #set par(first-line-indent: 0pt, justify: false)
    #align(left + horizon)[
$if(title)$
      #text(size: 24pt, weight: "bold")[$title$]
$endif$
$if(subtitle)$
      #v(6mm)
      #text(size: 16pt)[$subtitle$]
$endif$
$if(by-author)$
      #v(12mm)
      // Names boxed so the show rule hiding the template's author string does
      // not match, and so names do not break across lines
      #text(size: 12pt)[$for(by-author)$#box[$it.name.literal$]$sep$, $endfor$]
$endif$
$if(date)$
      #v(8mm)
      #text(size: 11pt)[$date$]
$endif$
    ]
  ],
$if(by-author)$
  author: "$for(by-author)$$it.name.literal$$sep$, $endfor$",
  copyright: [
    #set par(first-line-indent: 0pt, justify: false, spacing: 1em)
    #let equal-contributors = ($for(by-author)$$if(it.attributes.equal-contributor)$true, $endif$$endfor$)
    #text(weight: "bold")[$for(by-author)$#box[$it.name.literal$]#super[$for(it.affiliations)$$it.number$$sep$,$endfor$$if(it.attributes.equal-contributor)$#sym.dagger$endif$$if(it.attributes.corresponding)$\*$endif$]$sep$, $endfor$]

$for(by-affiliation)$
    #super[$it.number$] $if(it.department)$$it.department$, $endif$$it.name$$if(it.city)$, $it.city$$endif$$if(it.country)$, $it.country$$endif$ \
$endfor$

    #if equal-contributors.len() > 0 [#super[#sym.dagger] These authors contributed equally to this work. \ ]
$for(by-author)$
$if(it.attributes.corresponding)$
    #super[\*] Corresponding author: $it.name.literal$$if(it.email)$, $it.email$$endif$ \
$endif$
$endfor$
  ],
$endif$
$if(date)$
  date: "$date$",
$endif$
$if(lang)$
  lang: "$lang$",
$endif$
  main-color: accent-color,
  logo: {
    let logo-info = brand-logo.at("medium", default: none)
    if logo-info != none { image(logo-info.path, alt: logo-info.at("alt", default: none)) }
  },
$if(toc-depth)$
  outline-depth: $toc-depth$,
$endif$
$if(lof)$
  list-of-figure-title: "$if(crossref.lof-title)$$crossref.lof-title$$else$$crossref-lof-title$$endif$",
$endif$
$if(lot)$
  list-of-table-title: "$if(crossref.lot-title)$$crossref.lot-title$$else$$crossref-lot-title$$endif$",
$endif$
$if(margin-geometry)$
  padded-heading-number: false,
$endif$
  )(doc)
}

$if(margin-geometry)$
// Configure marginalia page geometry for book context
// Geometry computed by Quarto's meta.lua filter (typstGeometryFromPaperWidth)
// IMPORTANT: This must come AFTER book.with() to override the book format's margin settings
#import "@preview/marginalia:0.3.1" as marginalia

#show: marginalia.setup.with(
  inner: (
    far: $margin-geometry.inner.far$,
    width: $margin-geometry.inner.width$,
    sep: $margin-geometry.inner.separation$,
  ),
  outer: (
    far: $margin-geometry.outer.far$,
    width: $margin-geometry.outer.width$,
    sep: $margin-geometry.outer.separation$,
  ),
  top: $if(margin.top)$$margin.top$$else$1.25in$endif$,
  bottom: $if(margin.bottom)$$margin.bottom$$else$1.25in$endif$,
  // CRITICAL: Enable book mode for recto/verso awareness
  book: true,
  clearance: $margin-geometry.clearance$,
)
$endif$
