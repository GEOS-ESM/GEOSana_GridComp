file(GLOB perlscripts *.pl)
file(GLOB pythonscripts *.py)

install(
   FILES ${perscripts} ${pythonscripts}
   DESTINATION bin)
