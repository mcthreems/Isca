netcdf lgam40.q{                                                                

dimensions:
    lat            =   1;
    lon            =   1;
    plev           =  40;


variables:
    float lat(lat);                                                           
             lat:units = "degree";                                            
             lat:title = "LATITUDE";                                         
    float lon(lon);                                                           
             lon:units = "degree";                                            
             lon:title = "LONGITUDE";                                        
    float plev(plev);                                                         
             plev:units = "Pa";                                               
             plev:title = "PRESSURE";                                        

    float q(plev,lon,lat);                                                     
             q:units = "None";                                                
             q:title = "MMR OF WATER VAPOUR";                                                

data:                                                                           
              lat =   .000000E+00;
              lon =   .000000E+00;
             plev =   .252225E+03,  .754675E+03,  .125713E+04,  .175958E+04,
                      .226203E+04,  .276448E+04,  .326693E+04,  .376938E+04,
                      .427183E+04,  .477428E+04,  .527673E+04,  .577918E+04,
                      .628163E+04,  .678408E+04,  .728653E+04,  .778898E+04,
                      .829143E+04,  .879388E+04,  .929633E+04,  .979878E+04,
                      .122405E+05,  .166215E+05,  .210025E+05,  .253835E+05,
                      .297645E+05,  .341455E+05,  .385265E+05,  .429075E+05,
                      .472885E+05,  .516695E+05,  .560505E+05,  .604315E+05,
                      .648125E+05,  .691935E+05,  .735745E+05,  .779555E+05,
                      .823365E+05,  .867175E+05,  .910985E+05,  .954795E+05;
                q =   .328117E-05,  .311570E-05,  .304693E-05,  .300093E-05,
                      .296702E-05,  .294022E-05,  .287730E-05,  .279197E-05,
                      .271942E-05,  .265653E-05,  .259363E-05,  .253202E-05,
                      .247682E-05,  .242695E-05,  .239334E-05,  .237092E-05,
                      .235009E-05,  .233066E-05,  .231245E-05,  .229534E-05,
                      .329657E-05,  .721793E-05,  .194423E-04,  .513167E-04,
                      .174448E-03,  .315632E-03,  .518584E-03,  .757840E-03,
                      .103505E-02,  .136361E-02,  .173368E-02,  .216469E-02,
                      .266116E-02,  .322759E-02,  .394716E-02,  .479237E-02,
                      .575715E-02,  .687354E-02,  .817614E-02,  .964662E-02;

}                                                                               