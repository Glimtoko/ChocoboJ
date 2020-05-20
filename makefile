target = ChocoboCompiled/bin/Chocobo
installloc = /prod/ChocoboJ/bin
installlink = /prod/bin

# Always force a rebuild
.PHONY: $(target)

$(target):
	julia compile.jl

install: $(target)
	cp $(target) $(installloc)/chocoboj
	ln -sf $(installloc)/chocoboj $(installlink)

