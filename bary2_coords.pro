pro bary2_coords, msun, mass, x, y, z, vx, vy, vz, x_sun, y_sun, z_sun, $
  vx_sun, vy_sun, vz_sun, xb, yb, zb, vxb, vyb, vzb

;Convert heliocentric coordinates into barycentric coordinates. Inputs
;are the sun's mass and the planets' masses & coordinates. Outputs
;are the sun's and the planets' barycentric coordinates.
;
;Joe Hahn
;jhahn@spacescience.org

;Number of planets and time samplings
sz=size(x)
Np=sz(1)
Nt=sz(2)

;Sun's barycentric coordinates
x_sun=dblarr(Nt)
y_sun=dblarr(Nt)
z_sun=dblarr(Nt)
vx_sun=dblarr(Nt)
vy_sun=dblarr(Nt)
vz_sun=dblarr(Nt)
for tm=0, Nt-1 do begin $
  mt=total(mass(*,tm))+msun
  x_sun(tm)=-total(mass(*,tm)*x(*,tm))/mt & $
  y_sun(tm)=-total(mass(*,tm)*y(*,tm))/mt & $
  z_sun(tm)=-total(mass(*,tm)*z(*,tm))/mt & $
  vx_sun(tm)=-total(mass(*,tm)*vx(*,tm))/mt & $
  vy_sun(tm)=-total(mass(*,tm)*vy(*,tm))/mt & $
  vz_sun(tm)=-total(mass(*,tm)*vz(*,tm))/mt & $
endfor



;Planets' barycentric coordinates
xb=dblarr(Np, Nt)
yb=dblarr(Np, Nt)
zb=dblarr(Np, Nt)
vxb=dblarr(Np, Nt)
vyb=dblarr(Np, Nt)
vzb=dblarr(Np, Nt)
for tm=0, Nt-1 do begin $
  xb(0,tm)=x(*,tm)+x_sun(tm) & $
  yb(0,tm)=y(*,tm)+y_sun(tm) & $
  zb(0,tm)=z(*,tm)+z_sun(tm) & $
  vxb(0,tm)=vx(*,tm)+vx_sun(tm) & $
  vyb(0,tm)=vy(*,tm)+vy_sun(tm) & $
  vzb(0,tm)=vz(*,tm)+vz_sun(tm) & $
endfor

return
end
