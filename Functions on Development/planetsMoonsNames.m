function [namepl, namecentral] = planetsMoonsNames(idpl, idcentral)

if idcentral == 1 % --> SS planets
    namecentral = "Sun";
    names = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Uranus", "Neptune", "Pluto"];
elseif idcentral == 5 % --> Jupiter moons
    namecentral = "Jupiter";
    names = [ "Io", "Europa", "Ganymede", "Callisto" ];
elseif idcentral == 6 % --> Saturn moons
    namecentral = "Saturn";
    names = [ "Enceladus", "Tethys", "Dione", "Rhea", "Titan" ];
elseif idcentral == 7 % --> Uranus moons
    namecentral = "Uranus";
    names = [ "Miranda", "Ariel", "Umbriel", "Titania", "Oberon" ];
end

namepl = char(names(idpl));

end