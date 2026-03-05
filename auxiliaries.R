# Lists of groups
Stages.vec = c('Acute','Recovery','Convalescent', 'Extended')
Stages.HC.vec = c('Acute','Recovery','Convalescent', 'Extended',"Healthy Volunteers")

# Lists of comparisons ----

stageTest.paired = combn(Stages.vec, 2, simplify = FALSE)
stageTest.unpaired = list(
  c('Acute','Healthy Volunteers'),
  c('Recovery','Healthy Volunteers'),
  c('Convalescent','Healthy Volunteers'),
  c('Extended','Healthy Volunteers')
)

# Lists of color schemes -----

colorscheme <- c(
  "Acute"                   = "#b2182b",  # deep red
  "Recovery"                = "#ef8a62",  # reddish-orange
  "Convalescent"   = "#fdbb84",  # warm peach
  "Extended"   = "#7f7f7f",  # light tan (still visible)
  "Healthy Volunteers"       = "#000000"   # solid medium gray
)

colorscheme_omicron <- c(
  "Pre-Omicron" = pal_jco("default")(10)[1],  
  "Omicron" = pal_jco("default")(10)[2]  
)

colorscheme_severity <- c(
  "Hypoxic"  = pal_cosmic('hallmarks_light')(5)[2],
  "Normoxic" = pal_cosmic('hallmarks_light')(5)[5]
)

colorscheme_antibio <- c(
  "Antiviral + Antibiotics" = "orange",
  "Not used" = "darkblue",
  "Antiviral" = "darkred"
)

colorscheme_treatment <- c(
  "Treated" = '#C1440E',
  "Untreated" = "darkblue"
)


colorscheme_vaccination <- c(
  "Vaccinated" = '#1f77b4',
  "Unvaccinated" = "#ff7f0e"
)

colorscheme_proteome <- c(
  "Available"   = pal_cosmic('hallmarks_light')(5)[3],
  "Unavailable" = pal_cosmic('hallmarks_light')(5)[1]
)