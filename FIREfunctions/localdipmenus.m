function [menulist,excludelist] = localdipmenus(menulist)
I = size(menulist,1)+1;
menulist{I,1} = 'FIRE';
menulist{I,2} = {'fire_ims','fire_locs','-','frc','flc','fpc','-'};
excludelist = {'isect','sphere_tesselation','wkk','gfca3D_sub','frctoresolution'};