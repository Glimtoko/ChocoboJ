target = ChocoboCompiled/bin/Chocobo
installloc = /prod/ChocoboJ/bin
installlink = /prod/bin

# Always force a rebuild
.PHONY: install

# Source Files
source = Chocobo/src/Chocobo.jl \
         Chocobo/src/data.jl \
         Chocobo/src/eos.jl \
         Chocobo/src/mainloop.jl \
         Chocobo/src/mesh_reader.jl \
         Chocobo/src/new_geometry.jl \
         Chocobo/src/new_hydro.jl \
         Chocobo/src/old_geometry.jl \
         Chocobo/src/output_silo.jl \
         Chocobo/src/output_text.jl \
         Chocobo/src/set_initial_conditions.jl \
         Chocobo/src/user_input.jl \


$(target): $(source)
	julia compile.jl

install: $(target)
	-mkdir -p $(installloc)
	-mkdir -p $(installlink)
	cp $(target) $(installloc)/chocoboj
	ln -sf $(installloc)/chocoboj $(installlink)

