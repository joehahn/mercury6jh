pro angular_momentum, msun, mass, x_sun, y_sun, z_sun, vx_sun, vy_sun, vz_sun,$
  xb, yb, zb, vxb, vyb, vzb, L, Lx=Lx, Ly=Ly, Lz=Lz

;Calculate the system's total angular momentum about the barycenter.
;Inputs are the sun's mass, the planets' masses, the Sun's barycentric 
;coordinates, and planet's barycentric coordinates. Mass units are
;solar masses, length units are AU, and velocity units are 2*PI*AU/year,
;with the gravitational constant $G$ assumed to be unity. Output is the 
;system's total angular momentum, in units of solar masses*2*PI*AU^2/year.
;Optional outputs are the xyz components of the total angular momentum.

;Number of time samplings
sz=size(xb)
Nt=sz(2)

;Components of the planets' angular momenta
Lxi=mass*(yb*vzb-zb*vyb)
Lyi=mass*(zb*vxb-xb*vzb)
Lzi=mass*(xb*vyb-yb*vxb)

;Components of the Sun's angular momentum
Lx_sun=msun*(y_sun*vz_sun-z_sun*vy_sun)
Ly_sun=msun*(z_sun*vx_sun-x_sun*vz_sun)
Lz_sun=msun*(x_sun*vy_sun-y_sun*vx_sun)

;Components of the total angular momentum
Lx=dblarr(Nt)
Ly=dblarr(Nt)
Lz=dblarr(Nt)
for tm=0, Nt-1 do begin $
  Lx(tm)=total(Lxi(*,tm))+Lx_sun(tm) & $
  Ly(tm)=total(Lyi(*,tm))+Ly_sun(tm) & $
  Lz(tm)=total(Lzi(*,tm))+Lz_sun(tm) & $
endfor

;total angular momentum
L=sqrt(Lx^2+Ly^2+Lz^2)

return
end


