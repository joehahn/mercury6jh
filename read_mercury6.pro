;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Read the compressed Mercury6 output. 
;	filename=string containing input filename		;
;	t=time array in days					;
;	x,y,z=2D arrays containing the particles'		;
;		cartesian positions with indices		;	
;		x(particle, time) in units of AU.		;
;	vx,vy,vz=2D velocity arrays in units of AU/day		;
;	optional output name=string array containing each	;
;		particle's name.				;
;	set twopiyears=1 for t/x/vx in units of years, AU, and	;	
;		2*PI*AU/year					;
;	Np_max=maximum number of particles (default=1111l)	;
;	Nt_max=maximum number of time-outputs (default=1111l)	;
;								;
;To do: read the particles' spin data.				;
;								;
;Usage:								;
;								;
;read_mercury6,filename,t,Mcen,mass,x,y,z,vx,vy,vz, $		;
;  name=name,twopiyears=1,Np_max=1111l,Nt_max=1111l		;
;								;
;By: Joe Hahn, jhahn@spacescience.org				;
;								;
;Version 1.0							;
;February 28, 2008							;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;The following functions are literal translations of John Chamber's mio_c*	;
;functions that are found in his element6.f file that are also called 
;by read_mercury6, which is the final procedure in this file.		;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function mio_c2_i,c

N=2
c_byte=bytarr(N)
for j=0,N-1 do begin $
  c_byte(j)=byte(strmid(c,j,1)) & $
  if (c_byte(j) ge 26) then c_byte(j)=c_byte(j)-1 & $
  if (c_byte(j) ge 13) then c_byte(j)=c_byte(j)-4 & $
endfor
result =  1.5872763924382153d-5 * (251*c_byte(0) + c_byte(1))
return,result
end



function mio_c2re, c, xmin, xmax, nchar

y=0.d0
for j=nchar-1,0,-1 do begin $
  cj_byte=byte(strmid(c,j,1)) & $
  cj_byte=cj_byte(0) & $
  y=(y+double(cj_byte-32))/224.d0 & $
endfor
result=xmin+y*(xmax-xmin)
return,result
end



function mio_c2fl, c

x=mio_c2re(c,0.d0,1.d0,7)
x=x*2.d0-1.d0
ex=byte(strmid(c,7,1))
ex=ex(0)
ex=double(ex-32-112)
result=x*(10.d0^ex)
return,result
end



function mio_c4_i, c

N=4
c_byte=bytarr(N)
for j=0,N-1 do begin $
  c_byte(j)=byte(strmid(c,j,1)) & $
  if (c_byte(j) ge 26) then c_byte(j)=c_byte(j)-1 & $
  if (c_byte(j) ge 13) then c_byte(j)=c_byte(j)-4 & $
endfor
result = 2.5194463459916754d-10 * $
         (63001*c_byte(1) + 251*c_byte(2) + c_byte(3)) $
       + 3.9840637450199202d-3 * c_byte(0)
return,result
end



function mio_c7_i, c

N=7
c_byte=bytarr(N)
for j=0,N-1 do begin $
  c_byte(j)=byte(strmid(c,j,1)) & $
  if (c_byte(j) ge 26) then c_byte(j)=c_byte(j)-1 & $
  if (c_byte(j) ge 13) then c_byte(j)=c_byte(j)-4 & $
endfor
result = 1.59325d-17 * c_byte(6) $
  + 3.9990577071d-15 * $
    (63001*c_byte(3) + 251*c_byte(4) + c_byte(5)) $
  + 6.323810328439105d-8 * $
    (63001*c_byte(0) + 251*c_byte(1) + c_byte(2))
return,result
end



function mio_c8_r,c

c07=strmid(c,0,7)
c8=strmid(c,7,1)
x=mio_c7_i(c07)
x=x*2.d0-1.d0
ex=byte(c8)
ex=ex(0)
if (ex ge 26) then ex=ex-1
if (ex ge 13) then ex=ex-4
ex=ex-125
result=x*(10.0d^ex)
return,result
end



pro mco_ov2x,rcen,rmax,mcen,m,fr,theta,phi,fv,vtheta,vphi, $
  x,y,z,u,v,w

;Converts mercury6 output variables for an object to
;coordinates and velocities.
r = rcen * (10.d0^fr)
temp = sqrt(.5d0*(1.d0/fv - 1.d0))
v1 = sqrt(2.d0 * temp * (mcen + m) / r)
x = r * sin(theta) * cos(phi)
y = r * sin(theta) * sin(phi)
z = r * cos(theta)
u = v1 * sin(vtheta) * cos(vphi)
v = v1 * sin(vtheta) * sin(vphi)
w = v1 * cos(vtheta)
return
end



pro read_mercury6,filename,t,Mcen,mass,x,y,z,vx,vy,vz, $
  name=name,twopiyears=twopiyears,Np_max=Np_max,Nt_max=Nt_max

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Read the compressed Mercury6 output. Be sure to		;
;.run element6.pro before first call to this proceedure.	;
;	filename=string containing input filename		;
;	t=time array in days					;
;	x,y,z=2D arrays containing the particles'		;
;		cartesian positions with indices		;	
;		x(particle, time) in units of AU.		;
;	vx,vy,vz=2D velocity arrays in units of AU/day		;
;	optional output name=string array containing each	;
;		particle's name.				;
;	set twopiyears=1 for t/x/vx in units of years, AU, and	;	
;		2*PI*AU/year					;
;	Np_max=maximum number of particles (default=1111l)	;
;	Nt_max=maximum number of time-outputs (default=1111l)	;
;								;
;To do: read the particles' spin data.				;
;								;
;By: Joe Hahn, hahn@spacescience.org				;
;								;
;Version 1.0							;
;February 28, 2008							;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;anticipated maximum number of particles and outputs
if (keyword_set(Np_max) eq 0) then Np_max=1111l
if (keyword_set(Nt_max) eq 0) then Nt_max=1111l

;other constants
PI    = 3.141592653589793d0
TWOPI = 2.d0*PI
YEAR=365.25d0

;read input file and store in the string array data()
print,'Reading input file ',filename
data=strarr(2*Np_max*Nt_max)
get_lun,unit
openr,unit,filename
Nlines=0l
line=''
while not eof(unit) do begin $
  readf,unit,line & $
  data(Nlines)=line & $
  Nlines=Nlines+1l & $
endwhile
close,unit
free_lun,unit
data=data(0:Nlines-1)

;open file and read first line
line_no=0l
tm=0
print,'Decompressing data.'
while (line_no lt Nlines) do begin $

  ;read line and get data type and precision
  line=data(line_no) & $
  line_no=line_no+1l & $
  check=strmid(line,0,1) & $
  style=strmid(line,1,1) & $
  type=strmid(line,2,1) & $
  line_length=strlen(line) & $

  ;Read special input (eg, objects names, code, and masses, etc)
  if (type eq 'a') then begin $
    ;get time, nbig, nsml, Mcen, J's, rcen, rmax, rfac
    algor=fix(strmid(line,3,2)) & $
    precision=fix(strmid(line,67,2)) & $
    cc=strmid(line,5,62) & $
    time_str=strmid(cc,0,8) & $
    time=mio_c2fl(time_str) & $		;units of days
    nbig_str=strmid(cc,8,8) & $
    nbig=long(0.5d0+mio_c2re(nbig_str,0.d0, 11239424.d0, 3)) & $
    nsml_str=strmid(cc,11,8) & $
    nsml=long(0.5d0+mio_c2re(nsml_str,0.d0, 11239424.d0, 3)) & $
    Mcen_str=strmid(cc,14,8) & $
    Mcen=mio_c2fl(Mcen_str) & $
    J1_str=strmid(cc,22,8) & $
    J1=mio_c2fl(J1_str) & $
    J2_str=strmid(cc,30,8) & $
    J2=mio_c2fl(J2_str) & $
    J3_str=strmid(cc,38,8) & $
    J3=mio_c2fl(J3_str) & $
    rcen_str=strmid(cc,46,8) & $
    rcen=mio_c2fl(rcen_str) & $
    rmax_str=strmid(cc,54,8) & $
    rmax=mio_c2fl(rmax_str) & $
    rfac=alog10(rmax/rcen) & $

    ;output precision
    if (precision eq 1) then nchar=2 & $
    if (precision eq 2) then nchar=4 & $
    if (precision eq 3) then nchar=7 & $

    ;create storage on first entry
    Np=nbig+nsml & $
    if (tm eq 0) then begin $
      t=dblarr(Nt_max) & $
      Np_max=Np & $
      mass=dblarr(Np_max,Nt_max) & $
      x=dblarr(Np_max,Nt_max) & $
      y=dblarr(Np_max,Nt_max) & $
      z=dblarr(Np_max,Nt_max) & $
      vx=dblarr(Np_max,Nt_max) & $
      vy=dblarr(Np_max,Nt_max) & $
      vz=dblarr(Np_max,Nt_max) & $
      name=strarr(Np_max) & $
    endif & $

    ;get particle names, revised index, and mass
    index=indgen(Np, long=1) & $
    for p=0L,Np-1 do begin $
      line=data(line_no) & $
      line_no=line_no+1l & $
      c=strmid(line,0,8) & $
      nme=strmid(line,3,8) & $
      c=strmid(line,11,8) & $
      mss=mio_c2fl(c) & $
      if (tm eq 0) then name(p)=nme & $
      if (tm ne 0) then index(p)=where(nme eq name) & $  
      mass(index(p),tm:*)=mss & $    
    endfor & $

  endif & $

  ;otherwise read normal input (eg., coordinates, etc)
  if (type eq 'b') then begin $

    ;read time, nbig, nsml
    cc=strmid(line,3,14) & $
    time_str=strmid(cc,0,8) & $
    t(tm)=mio_c2fl(time_str) & $
    cc=strmid(line,3,14) & $
    time_str=strmid(cc,0,8) & $
    t(tm)=mio_c2fl(time_str) & $
    nbig_str=strmid(cc,8,8) & $
    nbig=fix(0.5d0+mio_c2re(nbig_str,0.d0, 11239424.d0, 3)) & $
    nsml_str=strmid(cc,11,8) & $
    nsml=fix(0.5d0+mio_c2re(nsml_str,0.d0, 11239424.d0, 3)) & $
    Np=nbig+nsml & $

    ;read each particle's coordinates
    for p=0,Np-1 do begin $
      ;get spherical coordinates
      line=data(line_no) & $
      line_no=line_no+1l & $
      c=strmid(line,0) & $
      str=strmid(c,0,8) & $
      code=fix(0.5d0+mio_c2re(str,0.d0, 11239424.d0, 3))-1 & $
      str=strmid(c,3,8) & $
      fr=mio_c2re(str, 0.d0, rfac, nchar) & $
      str=strmid(c,3+nchar,8) & $
      theta=mio_c2re(str, 0.d0, PI, nchar) & $
      str=strmid(c,3+2*nchar,8) & $
      phi=mio_c2re(str, 0.d0, TWOPI, nchar) & $
      str=strmid(c,3+3*nchar,8) & $
      fv=mio_c2re(str, 0.d0, 1.d0, nchar) & $
      str=strmid(c,3+4*nchar,8) & $
      vtheta=mio_c2re(str, 0.d0, PI, nchar) & $
      str=strmid(c,3+5*nchar,8) & $
      vphi=mio_c2re(str, 0.d0, TWOPI, nchar) & $
      idx=index(p) & $
      mss=mass(idx,tm) & $
      ;convert to cartesian coordinates having units of
      ;AU and TWOPI*AU/YEAR
      if ((fr gt 0.d0) and (fv gt 0.d0)) then begin $
        mco_ov2x,rcen,rmax,mcen,mss,fr,theta,phi,fv,vtheta,vphi,$
          x1,y1,z1,vx1,vy1,vz1 & $
        x(idx,tm)=x1 & $
        y(idx,tm)=y1 & $
        z(idx,tm)=z1 & $
        vx(idx,tm)=vx1 & $
        vy(idx,tm)=vy1 & $
        vz(idx,tm)=vz1 & $
      endif & $
    endfor & $

    ;increment time index tm
    tm=tm+1 & $

  endif & $

endwhile

;trim storage arrays
Nt=tm
t=t(0:Nt-1)
mass=mass(*,0:Nt-1)
x=x(*,0:Nt-1)
y=y(*,0:Nt-1)
z=z(*,0:Nt-1)
vx=vx(*,0:Nt-1)
vy=vy(*,0:Nt-1)
vz=vz(*,0:Nt-1)

;convert velocities to AU/day as needed
if (keyword_set(twopiyears) eq 0) then begin $
  factor=TWOPI/YEAR & $
  vx=vx*factor & $
  vy=vy*factor & $
  vz=vz*factor & $
endif

;convert time to years as needed
if (keyword_set(twopiyears) eq 1) then t=t/YEAR

return
end






