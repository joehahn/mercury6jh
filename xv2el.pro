pro xv2el,Msun,mass,x,y,z,vx,vy,vz,a,e,i,O,w,M, $
  mercury=mercury

;Convert Cartesian coordinates (,x,y,z,vx,vy,vz) into orbit
;elements (a,e,i,O,w,M). The primary mass is Msun and the 
;secondary masses are mass. All masses are in solar units,
;input distances are in units of AU, and velocities are in
;AU/(2*!pi*year). Output elements are in AU and radians.
;The GM of the system is Msun unless the mercury keyword
;is set for which GM=Msun+mass.
;
;Joe Hahn
;jhahn@spacescience.org

sz=size(x)
Ndims=sz(0)
arr_code=sz(Ndims+1)
if (arr_code eq 4) then TINY=2.0e-7
if (arr_code eq 5) then TINY=2.0d-15

;inclination i
hx=y*vz-z*vy
hy=z*vx-x*vz
hz=x*vy-y*vx
hxy=sqrt(hx^2+hy^2)
h2=hx^2+hy^2+hz^2
h=sqrt(h2)
i=make_array(size=sz)
r=sqrt(x^2+y^2+z^2)
j=where(r gt TINY)
if (j(0) ne -1) then i(j)=acos(hz(j)/h(j))

;sum primary and secondary masses unless mercury keyword is set
mu=Msun
if (keyword_set(mercury)) then mu=Msun+mass

;eccentricity e
ex=make_array(size=sz)
if (j(0) ne -1) then $
  ex(j)=(vy(j)*hz(j)-vz(j)*hy(j))/mu(j)-x(j)/r(j)
ey=make_array(size=sz)
if (j(0) ne -1) then $
  ey(j)=(vz(j)*hx(j)-vx(j)*hz(j))/mu(j)-y(j)/r(j)
ez=make_array(size=sz)
if (j(0) ne -1) then $
  ez(j)=(vx(j)*hy(j)-vy(j)*hx(j))/mu(j)-z(j)/r(j)
e=sqrt(ex^2+ey^2+ez^2)

;ascending node O
O=make_array(size=sz)
if (j(0) ne -1) then O(j)=atan(hx(j),-hy(j))

;argument of periapse w
ec=ex*cos(O)+ey*sin(O)
es=make_array(size=sz)
j=where(i gt TINY)
if (j(0) ne -1) then es(j)=ez(j)/sin(i(j))
w=atan(es,ec)

;semimajor axis a
v2=vx^2+vy^2+vz^2
a=make_array(size=sz)
j=where(r gt TINY)
if (j(0) ne -1) then a(j)=1.d0/(2.d0/r(j)-v2(j)/mu(j))

;mean anomaly M
M=make_array(size=sz)

;for elliptic orbits
rv=x*vx+y*vy+z*vz
k=where(rv lt 0.d0)
ec=make_array(size=sz)
ec(j)=1.d0-r(j)/a(j)
j=where(a gt 0.d0)
if (j(0) ne -1) then begin $ 
  es(j)=rv(j)/sqrt(abs(a(j))*mu(j)) & $
  Ea=make_array(size=sz) & $
  Ea(j)=atan(es(j),ec(j)) & $
  M(j)=Ea(j)-es(j) & $
endif

;for hyperbolic orbits
j=where(a lt 0.d0)
if (j(0) ne -1) then begin
  ZZ=make_array(size=sz)
  ZZ(j)=ec(j)/(e(j)>1.d0)
  F=make_array(size=sz)
  F(j)=alog(abs(ZZ(j)+sqrt(abs(ZZ(j)^2-1.d0))))
  if (k(0) ne -1) then F(k)=-F(k)
  M(j)=e(j)*sinh(F(j))-F(j)
endif

;for parabolic orbits
j=where((a eq 0.d0) and (r gt TINY))
if (j(0) ne -1) then begin
  q=make_array(size=sz)
  q(j)=0.5d0*h2(j)/mu(j)
  tau=make_array(size=sz)
  tau(j)=sqrt(r(j)/q(j)-1.d0)
  if (k(0) ne -1) then tau(k)=-tau(k)
  M(j)=sqrt(2.d0)*(tau(j)+(tau(j)^3)/3.d0)
endif

tp=6.283185307179586d0
j=where(O lt 0.d0)
if (j(0) ne -1) then O(j)=O(j)+tp
j=where(w lt 0.d0)
if (j(0) ne -1) then w(j)=w(j)+tp
j=where(M lt 0.d0)
if (j(0) ne -1) then M(j)=M(j)+tp

return
end
