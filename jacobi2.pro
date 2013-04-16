pro jacobi2, msun, x_sun, y_sun, z_sun, planet_idx, mass, xb, yb, zb, $
  vxb, vyb, vzb, J

;Calculates a massless particle's Jacobi integral. Inputs are the sun's 
;mass and barycentric coordinates, the planet's index , 
;and the masses and barycentric coordinates of all bodies. Mass units are
;solar masses, length units are AU, and velocity units are 2*PI*AU/year,
;with the gravitational constant $G$ assumed to be unity. Output is the 
;Jacobi integral for all bodies, in energy units, e.g.,
;solar masses*(2*PI*AU/year)^2.
;
;Joe Hahn
;jhahn@spacescience.org

;Number of bodies and time samplings
sz=size(xb)
Np=sz(1)
Nt=sz(2)

;planet's mass and coordinates
mp=mass(planet_idx, 0)
xp=xb(planet_idx, *)
yp=yb(planet_idx, *)
zp=zb(planet_idx, *)
vxp=vxb(planet_idx, *)
vyp=vyb(planet_idx, *)
vzp=vzb(planet_idx, *)
rp2=xp^2+yp^2+zp^2

;particles' helio- and planeto-centric distances
r_sun=dblarr(Np, Nt)
r_planet=dblarr(Np, Nt)
for tm=0, Nt-1 do begin $
  r_sun(0,tm)=sqrt( (xb(*,tm)-x_sun(tm))^2 + (yb(*,tm)-y_sun(tm))^2 + $
                    (zb(*,tm)-z_sun(tm))^2 ) & $
  r_planet(0,tm)=sqrt( (xb(*,tm)-xp(tm))^2 + (yb(*,tm)-yp(tm))^2 + $
                       (zb(*,tm)-zp(tm))^2 ) & $
endfor
r_planet(planet_idx, *)=1.d0

;particle's specific energy
v2=vxb^2+vyb^2+vzb^2
Energy=(0.5d0)*v2-msun/r_sun-mp/r_planet

;particle's specific angular momentum
Lx=yb*vzb-zb*vyb
Ly=zb*vxb-xb*vzb
Lz=xb*vyb-yb*vxb

;planet's specific angular momenta
Lxp=Lx(planet_idx, *)
Lyp=Ly(planet_idx, *)
Lzp=Lz(planet_idx, *)

;particle's Jacobi integral
J=Energy
for tm=0, Nt-1 do begin $
 J(0,tm)=J(*,tm) - $
   ( Lxp(tm)*Lx(*,tm) + Lyp(tm)*Ly(*,tm) + Lzp(tm)*Lz(*,tm))/rp2(tm) & $
endfor
J(planet_idx, *)=0.d0
J=J*(-2.d0)

return
end
