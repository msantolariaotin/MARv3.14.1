# MARv3.14.1

-------------------------------------------------------------------------------------------

#                                    MAR model                                          

-------------------------------------------------------------------------------------------

MAR (Modèle Atmosphérique Régional) is a regional climate model designed to simulate atmospheric processes and land–atmosphere interactions, with specific capabilities for applications over mountainous and polar regions.

Before compiling the MAR code, the integration domain, model grid and surface conditions must be defined to prepare the 6-hourly meteorological and surface forcing fields. This preprocessing step is performed using the nesting processor NESTOR (NESTing ORganization forcing) utility program (e.g., https://zenodo.org/records/15878925; see https://gitlab.uliege.be/groups/especes/mar/-/wikis/guide-domain-setup).

The MAR source code is available from https://gitlab.uliege.be/especes/mar/stable/, and the MAR source archives can be found here: http://ftp.climato.be/fettweis/MAR/.

Further documentation on the model and how to run it is available in the https://gitlab.uliege.be/groups/especes/mar/-/wikis/home/.

This README documents the specific modifications and configuration settings implemented in MAR and NESTOR for the simulations over High Mountain Asia (HMA) using MARv3.14.1. These modifications include, in particular, adaptations to the land-cover and ice mask, which are required to properly represent the regional characteristics of HMA. 

Two simuations have been already performed using this configuration at 7km x 7km horizontal resoluation during 1950-2025:

° The Central Himalaya (82.5°E–89.4°E, 26.6°N–30.5°N), https://zenodo.org/records/24011336

° The Pamir region (68.4°E–76.3°E, 36.8°N–40.7°N), https://zenodo.org/records/23921900

For version MARv3.14.0, please see https://zenodo.org/records/13151275.




-------------------------------------------------------------------------------------------

# Specific configurations for running MAR in the Northern Hemisphere over High Mountain Asia 
                                     
-------------------------------------------------------------------------------------------

* In INI.ctr :
 
WLDreg=GR for having permanent ice/snow (=GR, nsno=30,mw=2) instead of seasonal (=EU, nsno=10, mw=3). Having permanent ice/snow outside of ice sheets requires a mix between EU and GR with mw=2, ice in k=1 and dominant vegetation in k=2.


* In NST.ctr (activate land cover and vegetation) :


n22=T # land cover
n24=T #Globcov v2.2
n26=T #LAI/GLF from MERRA2 (50 km)

* In CODE.ctr :

LN=T #to activate LAI ("#define LN .true. /"*" Soil /Vegetation Model: LAI(x,y,t) prescribed(MARglf.DAT "*"/" >> $MAR_pp_def) for having LAI that varies with time

-----------------------------------
NESTOR settings
-----------------------------------

* In case of large domains (>500x200 grid points), in NESTOR/src/Compile : compile with -O3 -w -zero -mcmodel=large. Also consider using a high mzabso value (>= 8) and mz >= 30.

* In NESTOR/src/ : make sure that the dimensions of the NESTOR input files are correctly specified in NSTdim.inc (parameters ni,nj and nk) and adapt USRgrd.f to impose dominant vegetation instead of toundra in the subpixel and then allow subpixels to contain either ice, tundra or the dominant vegetation class.

* The section Topography for ETOPOg must contains :

'''
 IF (subnam.eq.'ETOPOg'.or.subnam.eq.'GLOcov') THEN

  write (*,*) 'Special ice mask for HMA Simulation'
  write (*,*)

  do i=1,mx ; do j=1,my

   !added from GLOcov
   if(NSTsvt(i,j,1)==-1.and.NSTsfr(i,j,1)>0) then
     NSTice(i,j)   = NSTsfr(i,j,1)
   else
     NSTice(i,j)   = 0.
   endif
   !end added from GLOcov

   if(NSTice(i,j)>10) then ! if ice > 10% => ice + tundra

    NSTsfr(i,j,1) = NSTice(i,j)
    NSTsfr(i,j,2) = min(100.,max(0.,100. - NSTsfr(i,j,1)))

    NSTsvt(i,j,1) =-1
    NSTsvt(i,j,2) = 1

   endif

   !if ice <= 10 % => vegetation only in k=2
   if(NSTice(i,j)<=10.and.NSTsvt(i,j,1)>=0) then ! only 1 vege type over Tundra

    if(i<=10.or.i>=mx-10.or.j<=10.or.j>=my-10) then
     if(NSTsvt(i,j,1)<NSTsvt(i,j,2)) NSTsvt(i,j,2)=NSTsvt(i,j,1)
    else
     if(NSTsfr(i,j,1)>NSTsfr(i,j,2)) NSTsvt(i,j,2)=NSTsvt(i,j,1)
    endif
    ! for having dominant vegetation  
    NSTsvt(i,j,1) = 1  ! tundra for k = 1 but NSTsfr(i,j,1) = 0 => only k=2

    NSTsfr(i,j,1) = 0
    NSTsfr(i,j,2) = min(100.,max(0.,100. - NSTsfr(i,j,1)))

   endif

  do k=1,mw
   NSTveg(i,j,k) = NSTsvt(i,j,k)
   NSTvfr(i,j,k) = NSTsfr(i,j,k)
  enddo

  enddo ; enddo

  ENDIF
'''


* Ice mask correction

We have also defined a new ice mask by correcting an overestimation in MAR/in/input/VEGE/glcesa3a.nc, which was caused by the misclassification of seasonal snow as glacier ice. To use this new file, the reading of input/ICEmask/ICEmask_full.nc in GLOcov.f must be commented out :

'''
! IF (region.eq."GRD".or.region.eq."ANT") THEN ! icemask = arrIce(kk,ll-minL) ! ENDIF
''''

and

'''
! IF ((region.eq."GRD".or.region.eq."ANT").and.icemask>0) THEN ! lc(-1)=lc(-1)+1.! Permanent snow and ice ! ELSE
'''

and the assignment of permanent snow and ice (cov==220) originally assigned to bare soil areas must be reassigned to permanent snow

'''
 !if(cov==220) lc(0)=lc(0)+1 ! Permanent snow and ice => Bare areas
    if(cov==220) lc(-1)=lc(-1)+1.
'''

* Vegetation correction

Similarly, vegetation class 11 Post-flooding or irrigated croplands (or aquatic) from the ESA GlobCover file glcesa3a.nc was mistakenly assigned to bare soils but must be reassigned to crops due to the inhibiting effect of dry bare soils on convection:

'''
 !if(cov==11) lc(0)=lc(0)+1 ! NO VEGETATION
    if(cov==11)  lc(3)=lc(3)+1.
'''

* The vertical grid in MARgrd.ctr has defined with :

'''

 2.                      | zzmin= STD NEW (0=>OK) Vertical discretisation 
 1.8                     | aavu=  STD NEW (0=>OK) "
 1.109                   | bbvu=  STD NEW (0=>OK) "
1000                    | ccvu=  STD NEW (0=>OK) "

'''
Note that with nk = 30 vertical levels and mzabso = 5 (defined in code/mar_dim_mod.f90), the last lowest level nudged in the free atmosphere lies at approximately 11 km asl.

* Humidity correction

There is a humidity correction in NESTOR (located in sim/GRq/input/NESTOR/src/NSTint.f) that is hard-coded and dependent on the horizontal resolution, as reported here :

'''

 C      WARNING increase of specific humidity
 C +    =====================================

                   correction = 1.005 ! by default +0.05 %

 if(NST_dx<=15000) correction = 1.01  ! correction minimum +1%
 if(NST_dx<=10000) correction = 1.05  ! correction minimum +5%
 if(NST_dx<=7500)  correction = 1.075 ! correction minimum +7.5%
 if(NST_dx<=5000)  correction = 1.10  ! correction minimum +10%

                   correction = min(1.2,max(0.8,correction))
''''

This correction is applied directly to the forcing files (outputs of NESTOR) and can be further adjusted in the code in marmagic_mod.f90. For example, at 5 km resolution, the correction by default is +10%.



-----------------------------------
MAR settings
-----------------------------------


*In code/sbcnew.f cut the reading of MARini files, use INI=4 in INIsnow and initialise the ice mask with either - 1) ice everywhere for sh < 7000 m with ice at the bottom and snow at the top of the snowpack: 

'''
if (INI==4)  	 .  rosSNo(i,j,n,k)   = min(920.,max(500., 500. + (7000.-sh(i,j)))) - 2) only ice at 920 kg/m3 if (INI==4)  	 .  rosSNo(i,j,n,k)   = 920.
'''

To cut the reading of MARini files and impose INI = 4 : if(TypeGL(1:2)=="GR") then if(mx==164 .and.my==116) TypeGL ='HMA' ! HMA

* In INIsnow if(TypeGL(1:2) == "HM") INI = 4

* Call INIsnow if(TypeGL(1:2)=="HM".or.TypeGL(1:2)=="AN") then IF(itexpe == 0) call INIsnow(TypeGL)

Check in the output files ICE*.nc after 1 day if there is 20 m of snow (ZN3[k=1,l=1]), vegetation (VEG[k=1,l=1]) et snow at the right density RO1[k=1,l=1] at 300 kg/m3 at the top and 920 kg/m3 at the bottom. 
Also check the surface type. For instance, to verify that  no surface lake (SRF=1) has been assigned in the middle of nowhere. If lake, check in lac ETOPO.f: at the end of the main routine, NST_sol is assigned type 1:sea if topography is <0.1 m


-----------------------------------
Recommendations
-----------------------------------

*For the values of humidity_magic, cloud_magic, and rain_snow_limit over HMA from X. Fettweis :

-humidity_magic = 10. The larger this value, the greater the diffusion, which results in more clouds and warmer conditions. A value of 10 is considered the median.

- cloud_magic = 0. The larger it is, the more rain/snow precipitation is converted back to cloud rain/ice, leads to more clouds and delayed precipitation. To bet set between 0 and 0.2 max.

-rain_snow_limit = 273.15. Since we are dealing with a tropical climate with significant convection, this is set at 273.15 K.


*mzabso can be increased to 8 or 9 for large domains to prevent MAR from drifting too much (as long as the highest pixel is outside the nudged area). If mz = 24, the first level should probably be set at 3m instead of 2m.


*mzhyd should remain as mzhyd = min(4, mzabso - 1) even for a large domain and high mzabso values


*cnos in hydmic.f90 is still a tuning parameter for adjusting the spatial distribution of precipitation, but probably not at first order. The most important thing is to have observations for proper validation. However, disregard observations in Nepal (behind the ridges) because the barrier effect will inevitably be underestimated:  in that case a higher cnos value is likely required.
