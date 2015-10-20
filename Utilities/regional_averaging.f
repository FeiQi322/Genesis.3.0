c       Calculate a regional spatial average (valr) of arr 
c       over the rectangle defined by alata, alatb, alona, and alonb.
c       Enter latitudes alata,alatb as degrees N (-90 to 90), and 
c       longitudes alona,alonb as deg E from Greenwich (-180 to 180).
c       Nb: This code doesn't work if the longitudinal region spans
c       the dateline (180 E).

      alona = -10.
      alonb =  40.
      alata =  20..
      alata =  60.

      zlona = mod (alona + 360., 360.)
      if (zlona.gt.180.) zlona = zlona - 360.   
      zlonb = mod (alonb + 360., 360.)
      if (zlonb.gt.180.) zlonb = zlonb - 360.   

      valr = 0.
      vala = 0.
      do j=1,nlat
        if (alat(j).ge.alata .and. alat(j).le.alatb) then  
          do i=1,nlon
            zlon = mod (alon(i) + 360., 360.)
            if (zlon.gt.180.) zlon = zlon - 360.   
            if (zlon.ge.zlona .and. zlon.le.zlonb) then
              valr = valr + weilat(j)*arr(i,j)
              vala = vala + weilat(j)           
            endif
          enddo
        endif
      enddo
      if (vala.gt.0.) valr = valr/vala

