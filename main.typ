#import "@local/templates:1.0.0" : *

#show: template

// Document
#make-title("Linear Regression")

#block(width:100%)[
  #set quote(block: true)
  #quote(attribution: [Albert Einstein])["Everything should be made as simple as possible, but not simpler."]
]

// Main contents of the document
#show outline: it => {
  show heading.where(level: 1): it_in => [
    #set block(above: 1.4em, below: 1em)
    #set text(font: "Arial")
    #it_in.body \
  ]
  it
}
#show outline.entry.where(level:1): it => {
  v(12pt, weak:true)
  strong(it)
}
#outline(title:[Contents], indent: auto)


// Main contents
#include "linear-regression.typ"
#include "broken-assumptions/weak-exogeneity.typ"
#include "broken-assumptions/multicollinearity.typ"

// References
#bibliography("ref.bib", title:"References")