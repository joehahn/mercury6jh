;demo.pro

;demonstrate how to use various IDL routines to
;examine mercury6.2jh output
;
;Joe Hahn
;jhahn@spacescience.org


;Read mercury6.2jh output, which are heliocentric coordinates
read_mercury6, 'xv.out', t, msun, mass, x, y, z, $
  vx, vy, vz, twopiyears=1, name=name

;get the number of particles Np and time samplings Nt
sz=size(x)
Np=sz(1)
Nt=sz(2)

;Convert coordinates to elements. Note that lost particle's
;have their coordinates zeroed, so for convenience, put these
;particles on wide circular orbits before converting to elements
rmin=0.1
rmax=9999.9
r=sqrt(x^2+y^2+z^2)
pe=where(r lt rmin)
if (pe(0) ne -1) then begin $
  x(pe)=rmax & $
  vy(pe)=1.0/sqrt(x(pe)) & $
endif 
xv2el,msun,mass,x,y,z,vx,vy,vz,a,e,i,O,w,M,mercury=1

;display all particles' e versus a for all times
window,xs=550,ys=550,retain=2
plot,a,e,psym=3,xrange=[0,2.5],yrange=[0,0.7], ystyle=1

;Get system's barycentric coordinates
bary2_coords, msun, mass, x, y, z, vx, vy, vz, x_sun, y_sun, z_sun, $
  vx_sun, vy_sun, vz_sun, xb, yb, zb, vxb, vyb, vzb

;extract the planet's barycentric coordinates
planet_idx=0
mp=mass(planet_idx,*)
xp=xb(planet_idx,*)
yp=yb(planet_idx,*)
zp=zb(planet_idx,*)
vxp=vxb(planet_idx,*)
vyp=vyb(planet_idx,*)
vzp=vzb(planet_idx,*)

;Calculate the system's energy in the barycentric frame,
;as well as the fractional energy error dE
energy2, msun, mp, x_sun, y_sun, z_sun, vx_sun, vy_sun, $
  vz_sun, xp, yp, zp, vxp, vyp, vzp, energy
dE=abs(energy/energy(0)-1.d0)

;Angular momentum in the barycentric frame,
;as well as the fractional error dL
angular_momentum, msun, mp, x_sun, y_sun, z_sun, $
    vx_sun, vy_sun, vz_sun, $
    xp, yp, zp, vxp, vyp, vzp, L, Lx=Lx, Ly=Ly, Lz=Lz
dL=abs(L/L(0)-1.d0)

;Calculate all particles' Jacobi integrals
jacobi2, msun, x_sun, y_sun, z_sun, planet_idx, mass, xb, yb, zb, $
  vxb, vyb, vzb, J

;check energy and angular momentum conservation
plot,t,dE,ylog=1,yrange=[1.0e-14, 1.0e-8],ystyle=1, $
  xtitle='time t    (years)', ytitle='!7D!3E/E and !7D!3L/L', $
  charsize=1.3
oplot,t, dL,linestyle=2

;check Jacobi conservation for all particles
dJ_max=0.0
for p=1, Np-1 do begin $
  dJ=abs(J(p,*)/J(p,0)-1.d0) & $
  ;make sure particle wasn't removed from the system
  tm=where((r(p,*) gt rmin) and (r(p,*) lt rmax)) & $
  if (p eq 1) then plot,t(tm),dJ(tm),ylog=1,  $
    yrange=[1.0e-12, 1.0e-6],psym=3,xtitle='time t    (years)', $
    ytitle='!7D!3J/J',charsize=1.3 & $
  oplot,t(tm),dJ(tm),psym=3 & $
  if (max(dJ(tm)) gt dJ_max) then dJ_max=max(dJ(tm)) & $
endfor
print,'J conserved to ', dJ_max
