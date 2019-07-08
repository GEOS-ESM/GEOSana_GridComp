file(GLOB perlscripts *.pl)
file(GLOB pythonscripts *.py)

install(
   PROGRAMS ${perscripts} ${pythonscripts}
   DESTINATION bin)
