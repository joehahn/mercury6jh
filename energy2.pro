pro energy2, msun, mass, x_sun, y_sun, z_sun, vx_sun, vy_sun, vz_sun, $
  xb, yb, zb, vxb, vyb, vzb, energy

;Calculate the system's total energy in the barycentric reference frame.
;Inputs are the sun's mass, the planets' masses, the Sun's barycentric 
;coordinates, and planet's barycentric coordinates. Mass units are
;solar masses, length units are AU, and velocity units are 2*PI*AU/year,
;with the gravitational constant $G$ assumed to be unity. Output is the 
;system's total energy, in units of solar masses*(2*PI*AU/year)^2. 
;
;Joe Hahn
;jhahn@spacescience.org

;Number of planets and time samplings
sz=size(xb)
Np=sz(1)
Nt=sz(2)

;Kinetic energy
KE=dblarr(Nt)
for tm=0, Nt-1 do begin $
  ;KE due to planets
  v2=vxb(*,tm)^2+vyb(*,tm)^2+vzb(*,tm)^2 & $
  KE(tm)=(0.5d0)*total(mass(*,tm)*v2) & $
  ;add KE due to Sun
  v2_sun=vx_sun(tm)^2+vy_sun(tm)^2+vz_sun(tm)^2 & $
  KE(tm)=KE(tm)+(0.5d0)*msun*v2_sun & $
endfor

;Will need the transpose of the Sun's coordinates
xt_sun=transpose(x_sun)
yt_sun=transpose(y_sun)
zt_sun=transpose(z_sun)

;Potential energy
PE=dblarr(1, Nt)
for j=0, Np-1 do begin $
  for k=j+1, Np-1 do begin $
    ;contributions from planet-planet interactions
    delta=sqrt( (xb(j,*)-xb(k,*))^2 + (yb(j,*)-yb(k,*))^2 + $
                (zb(j,*)-zb(k,*))^2 )
    Ujk=-mass(j,*)*mass(k,*)/delta & $
    PE=PE+Ujk & $
  endfor
  ;sun-planet interactions
  delta=sqrt( (xb(j,*)-xt_sun)^2 + (yb(j,*)-yt_sun)^2 + $
              (zb(j,*)-zt_sun)^2 ) & $ 
  Ujsun=-mass(j,*)*msun/delta & $
  PE=PE+Ujsun & $
endfor
PE=reform(PE)

;Total energy
energy=KE+PE

return
end
