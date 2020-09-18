within ;
package boiler
  package streams
    connector dryair
      "connector for connecting dry air streams"

       Real f(unit="kg/s") "Mass flow rate";

      Real P(
        unit="N/m2") "Stream Pressure";

      Real T(
        unit="K") "Temperature";




       annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Polygon(
              points={{-100,92},{-100,-92},{96,0},{-100,92}},
              lineColor={28,108,200},
              fillColor={238,46,47},
              fillPattern=FillPattern.Solid)}), Diagram(coordinateSystem(
              preserveAspectRatio=false), graphics={Polygon(
              points={{-100,92},{-100,-92},{96,0},{-100,92}},
              lineColor={28,108,200},
              fillColor={238,46,47},
              fillPattern=FillPattern.Solid)}));
    end dryair;

    connector water
       "connector for connecting water streams"
       Real f(
        unit="kg/s") "Mass flow rate";
      Real z_water(
        unit="") "mass fraction of water";
      Real z_steam(
        unit="")  "mass fraction of steam";
      Real P(
        unit="N/m2") "Stream Pressure";
      Real T(
        unit="K") "Temperature of the stream";
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Polygon(
              points={{-100,98},{-100,-98},{96,0},{-100,98}},
              lineColor={28,108,200},
              fillColor={28,108,200},
              fillPattern=FillPattern.Solid)}), Diagram(coordinateSystem(
              preserveAspectRatio=false), graphics={Polygon(
              points={{-100,98},{-100,-98},{96,0},{-100,98}},
              lineColor={28,108,200},
              fillColor={28,108,200},
              fillPattern=FillPattern.Solid)}));
    end water;
  end streams;

  package streamsources

    model watersource

      streams.water water
        annotation (Placement(transformation(extent={{-98,-102},{106,102}})));
            parameter Real Water=1 "fix value of Water in source stream";
      parameter Real Steam=0 "fix value of Steam in source stream";
      parameter Real T=30+273 "fix value of T in source stream";
       parameter Real P=1e5 "fix value of P in source stream";
       parameter Real f=100 "fix";


    equation

      water.z_water = Water "Water source eq";
      water.z_steam = Steam "Steam source eq";
      water.T = T "Temp source eq";
      water.f=f "f source eqn";
      water.P = P "P source eq";

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end watersource;

    model dryairsource
      streams.dryair dryair
        annotation (Placement(transformation(extent={{-100,-106},{102,106}})));

     parameter Real P=1e5 "Fix value of P in source stream";
      parameter Real T=300+273 "Fix value of T in source stream";
      parameter Real f=100 "flow rate";

    equation

      dryair.P = P "P source eq";
    // dryair.T = T "T source eq";
      if (time<5000) then
        dryair.T = T "T source eq";

      else
         dryair.T = T -20 "T source eq";
      end if;

    dryair.f=f "f source eqn";

      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end dryairsource;
  end streamsources;

  package equipments
    block economizer
    ///////////////////////inputs////
    input  streams.dryair dryair_in
        annotation (Placement(transformation(extent={{-70,36},{-50,56}})));

     input streams.water water_in "ambient water inlet to the economizer"
        annotation (Placement(transformation(extent={{64,0},{44,20}})));

        ////////////////////outputs/////
     output streams.water water_out "outlet water from the economizer"
        annotation (Placement(transformation(extent={{-56,0},{-76,20}})));
      output   streams.dryair dryair_out
        annotation (Placement(transformation(extent={{46,38},{66,58}})));
                          ///////////////output variables/////////////

    output Real Q(unit="J/s");
    output Real C_H(start=100000)
                                 "Mass Heat Capacity of hot stream ";
    output Real C_C(start=100000)
                                 "Mass Heat Capacity of cold stream ";
    output Real NTU;
    output Real epsilon "Efficiency of heat transfer";
    output Real Qmax(unit="J")
                              "Theoretical maximum heat transfer that can occur between the given streams";
    output Real ratio "Ratio of mass heat capacities of the streams exchanging heat";
    output Real C_min(start=500000, unit="J/K")  "Minimum of mass heat capacity of the streams exchanging heat";
    output Real C_max(start=5000000, unit="J/K") "Maximum of mass heat capacity of the streams exchanging heat";
    //output Real t1(unit="K");
    //output Real t2(unit="K");
    //output Real dt(unit="K");
                    ///////////////parameters/////////////////////

    parameter Real A(unit="m2") = 50 "Area of ht tr. m2";
    parameter Real rho_dryair(unit="kg/m3") = 1.2 "Density of dry air";
    parameter Real gamma(unit="m4") = 0.1  "fraction coefficient for flue gas flow";
    parameter Real fs(unit="m-4") = 2615 "friction coeff for steam flow";
    parameter Real m(unit="kg") = 100 "Mass Hold up";
    parameter Real U(unit="(J/s)/(m2.K)") = 25;
    parameter Real cp_dryair(unit="J/(kg.K)") = 1000;
    parameter Real cp_water(unit="J/(kg.K)") = 4180;

                  //////////////////equations//////////
    equation

                  //////////for dry air///////////

                  dryair_in.f = dryair_out.f  "Mass balance on flue gas flow";
    dryair_out.P = dryair_in.P  "dry air P loss eq";
     der(m*cp_dryair*dryair_out.T) = (dryair_in.f*cp_dryair*dryair_in.T -
        dryair_out.f*cp_dryair*dryair_out.T) - Q;

                  ////////////for Q/////////////////

      dryair_in.f*cp_dryair= C_H "Flue gas mass heat capacity";
      water_in.f*cp_water= C_C "Water gas mass heat capacity";
       if (C_H > C_C) then
        C_min = C_C "Min mass heat capacity";
        C_max = C_H "Max mass heat capacity";
      else
        C_min = C_H "Min mass heat capacity";
        C_max = C_C "Max mass heat capacity";
      end if;
      C_max*ratio = C_min "Heat capacity ratio";
      NTU*C_min = U*A "No of tr units";
      Qmax = C_min*(dryair_in.T - water_in.T)  "Max Heat tr. eq";
      Q = (epsilon*Qmax) "Actual heat tr";
      (1.0 - ratio*exp(-NTU*(1 - ratio)))*epsilon = (1.0 - exp(-NTU*(1 -
          ratio)))
                  "Heat tr efficiency eq";
    //t1=dryair_out.T-water_in.T;
    //t2=dryair_in.T-water_out.T;
    //dt=(t1-t2)/ln(t1/t2);
     //Q=U*A*dt;

                    /////////////for water////////////

     water_out.P = water_in.P "Water pressure";


      water_out.f = water_in.f "Mass balance on water";
      water_out.z_water = 1 "component mass balance";
      water_out.z_steam = 0 "component mass balance";

      der(m*cp_water*water_out.T) = (water_in.f*cp_water*water_in.T -
        water_out.f*cp_water*water_out.T) + Q;
    initial equation
        dryair_out.T=dryair_in.T;
       // dryair_out.T=300+273;
        //water_out.T=30+273;

     water_out.T=water_in.T;
             annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{-46,58},{40,0}},
              lineColor={28,108,200},
              fillColor={102,44,145},
              fillPattern=FillPattern.Solid)}),                             Diagram(
            coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
              extent={{-48,60},{42,0}},
              lineColor={28,108,200},
              fillColor={102,44,145},
              fillPattern=FillPattern.Solid)}));
    end economizer;

    block superheater

      ////////////inputs//////////////

     input streams.water steam_in
        annotation (Placement(transformation(extent={{80,-10},{60,10}})));
       input   streams.water attemperator_water
        annotation (Placement(transformation(extent={{10,-10},{-10,10}},
            rotation=-90,
            origin={44,-30})));
      input   streams.dryair dryair_in
        annotation (Placement(transformation(extent={{-80,52},{-60,72}})));

        ////////////outputs///////////

     output streams.water superheatedsteam_out
        annotation (Placement(transformation(extent={{-60,-10},{-80,10}})));

     output streams.dryair dryair_out
        annotation (Placement(transformation(extent={{60,50},{80,70}})));

     ///////////////output variables/////////////

    output Real Q(unit="J/s");
    output Real C_H(start=100000)
                                 "Mass Heat Capacity of hot stream ";
    output Real C_C(start=100000)
                                 "Mass Heat Capacity of cold stream ";
    output Real NTU;
    output Real epsilon "Efficiency of heat transfer";
    output Real Qmax(unit="J")
                              "Theoretical maximum heat transfer that can occur between the given streams";
    output Real ratio "Ratio of mass heat capacities of the streams exchanging heat";
    output Real C_min(start=500000, unit="J/K")  "Minimum of mass heat capacity of the streams exchanging heat";
    output Real C_max(start=5000000, unit="J/K") "Maximum of mass heat capacity of the streams exchanging heat";

                    ///////////////parameters/////////////////////

    parameter Real A(unit="m2") = 900 "Area of ht tr. m2";
    parameter Real rho_dryair(unit="kg/m3") = 1.2 "Density of dry air";
    parameter Real gamma(unit="m4") = 0.1  "fraction coefficient for flue gas flow";
    parameter Real fs(unit="m-4") = 2615 "friction coeff for steam flow";
    parameter Real m(unit="kg") = 100 "Mass Hold up";
    parameter Real U(unit="(J/s)/(m2.K)") = 25;
    parameter Real cp_dryair(unit="J/(kg.K)") = 1000;
    parameter Real cp_water(unit="J/(kg.K)") = 4180;
    parameter Real cp_steam(unit="J/(kg.K)") = 2000 "Specific Heat of Steam";


     //////////////////equations//////////
    equation

                  //////////for dry air///////////

    dryair_in.f = dryair_out.f  "Mass balance on flue gas flow";
    dryair_out.P = dryair_in.P  "dry air P loss eq";

     der(m*cp_dryair*dryair_out.T) =(dryair_in.f*cp_dryair*dryair_in.T - dryair_out.f*cp_dryair*
        dryair_out.T) - Q;



     ////////////for Q/////////////////

      dryair_in.f*cp_dryair= C_H "Flue gas mass heat capacity";
      steam_in.f*cp_water= C_C "Water gas mass heat capacity";
       if (C_H > C_C) then
        C_min = C_C "Min mass heat capacity";
        C_max = C_H "Max mass heat capacity";
      else
        C_min = C_H "Min mass heat capacity";
        C_max = C_C "Max mass heat capacity";
      end if;
      C_max*ratio = C_min "Heat capacity ratio";
      NTU*C_min = U*A "No of tr units";
      Qmax = C_min*(dryair_in.T - steam_in.T)  "Max Heat tr. eq";
      Q = (epsilon*Qmax) "Actual heat tr";
      (1.0 - ratio*exp(-NTU*(1 - ratio)))*epsilon = (1.0 - exp(-NTU*(1 -
          ratio)))
                  "Heat tr efficiency eq";

     /////////////for water////////////

     superheatedsteam_out.P = steam_in.P "Water pressure";


      superheatedsteam_out.f = steam_in.f + attemperator_water.f "Mass balance on water";
      superheatedsteam_out.z_water = 0 "component mass balance";
      superheatedsteam_out.z_steam = 1 "component mass balance";

      der(m*cp_steam*superheatedsteam_out.T) = (steam_in.f*cp_steam*steam_in.T +attemperator_water.f*cp_water*attemperator_water.T -
        superheatedsteam_out.f*cp_steam*superheatedsteam_out.T) + Q;

    initial equation
        dryair_out.T=dryair_in.T;
         superheatedsteam_out.T=steam_in.T;
    // dryair_out.T=600+273;
      //   superheatedsteam_out.T=100+273;





        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{-60,76},{60,-20}},
              lineColor={28,108,200},
              fillColor={244,125,35},
              fillPattern=FillPattern.Solid)}),                        Diagram(
            coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
              extent={{-60,76},{60,-20}},
              lineColor={28,108,200},
              fillColor={244,125,35},
              fillPattern=FillPattern.Solid)}));
    end superheater;

    block riser

      ////////////////inputs/////////////////

     input streams.dryair dryair_in annotation (Placement(transformation(extent={{-72,50},{-52,70}})));
     input streams.water water_in annotation (Placement(transformation(extent={{60,-10},{40,10}})));

      ///////////////////outputs//////////////

     output streams.dryair dryair_out annotation (Placement(transformation(extent={{40,50},{60,70}})));
     output streams.water water_out annotation (Placement(transformation(extent={{-50,-10},{-70,10}})));

    ///////////////output variables/////////////

    output Real Q(unit="J/s");
    output Real C_H(start=100000)
                                 "Mass Heat Capacity of hot stream ";
    output Real C_C(start=100000)
                                 "Mass Heat Capacity of cold stream ";
    output Real NTU;
    output Real epsilon "Efficiency of heat transfer";
    output Real Qmax(unit="J")
                              "Theoretical maximum heat transfer that can occur between the given streams";
    output Real ratio "Ratio of mass heat capacities of the streams exchanging heat";
    output Real C_min(start=500000, unit="J/K")  "Minimum of mass heat capacity of the streams exchanging heat";
    output Real C_max(start=5000000, unit="J/K") "Maximum of mass heat capacity of the streams exchanging heat";
     //output Real rho_riser(unit="kg/m3",start=560)  "Liq vapor mixture density at the riser";
          output Real sg(unit="kg/s");
      output Real rho_waterout(unit="kg/m3");
      output Real Hout(unit="J/kg");

                    ///////////////parameters/////////////////////

    parameter Real A(unit="m2") = 1000 "Area of ht tr. m2";
    parameter Real rho_dryair(unit="kg/m3") = 1.2 "Density of dry air";

    parameter Real rho_v(unit="kg/m3") = 4 "Density of water vapour ";
    parameter Real rho_w(unit="kg/m3") = 1000 "Density of liquid water";
    parameter Real gamma(unit="m4") = 0.1  "fraction coefficient for flue gas flow";
    parameter Real fs(unit="m-4") = 2615 "friction coeff for steam flow";
    parameter Real m(unit="kg") = 100
     "Mass Hold up";
    parameter Real U(unit="(J/s)/(m2.K)") = 25;
    parameter Real cp_dryair(unit="J/(kg.K)") = 1000;
    parameter Real cp_water(unit="J/(kg.K)") = 4180;
    parameter Real cp_steam(unit="J/(kg.K)") = 2000 "Specific Heat of Steam";
    parameter Real v_riser(unit="m3") = 7 "volume of ht tr. m3";
       parameter Real Evap_energy(unit="J/kg") = 2230000
        "Latent heat of vaporization";

     //////////////////equations//////////
    equation

                  //////////for dry air///////////

    dryair_in.f = dryair_out.f  "Mass balance on flue gas flow";
    dryair_out.P = dryair_in.P  "dry air P loss eq";
    der(m*cp_dryair*dryair_out.T) =(dryair_in.f*cp_dryair*dryair_in.T - dryair_out.f*cp_dryair*
        dryair_out.T) - Q;

    ////////////for Q/////////////////

      dryair_in.f*cp_dryair= C_H "Flue gas mass heat capacity";
      water_in.f*cp_water= C_C "Water gas mass heat capacity";
       if (C_H > C_C) then
        C_min = C_C "Min mass heat capacity";
        C_max = C_H "Max mass heat capacity";
      else
        C_min = C_H "Min mass heat capacity";
        C_max = C_C "Max mass heat capacity";
      end if;
      C_max*ratio = C_min "Heat capacity ratio";
      NTU*C_min = U*A "No of tr units";
      Qmax = C_min*(dryair_in.T - water_in.T)  "Max Heat tr. eq";
      Q = (epsilon*Qmax) "Actual heat tr";
      (1.0 - ratio*exp(-NTU*(1 - ratio)))*epsilon = (1.0 - exp(-NTU*(1 -
          ratio)))
                  "Heat tr efficiency eq";

     /////////////for water////////////
    water_out.P = water_in.P "P eq";
      rho_waterout = water_out.z_steam*rho_v + water_out.z_water*rho_w;

       water_in.f - water_out.f = v_riser*der(rho_waterout) "Mass balance on water";

    water_out.z_steam + water_out.z_water = 1;

    water_in.f*water_in.z_water - water_out.f*water_out.z_water - sg =der(v_riser*rho_waterout*water_out.z_water);
    //water_out.f = water_in.f "f eq new";

        Q + water_in.f*cp_water*water_in.T = Hout*water_out.f
        + der(v_riser*rho_waterout*Hout) "Heat balance eq";

        Hout=(water_out.f*cp_water*water_out.T + sg*2260000)/water_out.f;

       sg = (Q-water_out.f*cp_water*(water_out.T-water_in.T))/2260000;
      // sg = (Q-water_out.f*cp_water*(373-water_in.T)-water_out.f*cp_steam*(water_out.T-373))/2260000;
    //sg=1;
    //water_out.z_steam=0.2;

      //water_in.f*water_in.z_steam - water_out.f*water_out.z_steam + sg =v_riser*der(rho_v*water_out.z_steam);

    //water_in.f*water_in.z_water - water_out.f*water_out.z_water - sg =v_riser*der(rho_w*water_out.z_water);

    initial equation

     dryair_out.T=dryair_in.T;
     //water_out.z_water=1;
     //water_out.f=100;
     water_out.z_steam=0;
      water_out.T=water_in.T;
     // water_out.f=100;
    // dryair_out.T=500+273;
      // water_out.T=87+273;
     annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
              extent={{-48,76},{40,-24}},
              lineColor={28,108,200},
              fillColor={162,29,33},
              fillPattern=FillPattern.Solid)}),                      Diagram(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{-48,76},{40,-24}},
              lineColor={28,108,200},
              fillColor={162,29,33},
              fillPattern=FillPattern.Solid)}));
    end riser;

    block steamdrum

      //////////inputs///////////

     input streams.water water_in "water from econmizer"
        annotation (Placement(transformation(extent={{72,50},{52,70}})));
     input streams.water riserwater_in "water steam mixture from riser"
        annotation (Placement(transformation(extent={{72,-6},{52,14}})));

      ///////////outputs/////////

     output streams.water steam_out
        annotation (Placement(transformation(extent={{-44,50},{-64,70}})));
     output streams.water dwater_out
        annotation (Placement(transformation(extent={{-44,14},{-64,-6}})));

              output Real v_l(unit="m3");
              output Real v_v( unit="m3");


              parameter Real cp_water(unit="J/(kg.K)") = 4180;
    parameter Real cp_steam(unit="J/(kg.K)") = 2000;
              parameter Real V( unit = "m3")=10;
              parameter Real rho_w( unit= "kg/m3") = 1000;
              parameter Real rho_v( unit= "kg/m3") = 4;
              parameter Real f(unit="kg/s")= 5;
    // parameter Real F(unit="kg/s")= 6;
    equation

      steam_out.z_steam=1;
      steam_out.z_water=0;
      dwater_out.z_water=1;
      dwater_out.z_steam=0;
      steam_out.P=dwater_out.P;
      steam_out.P=water_in.P;

      ////////////MASS BALANCES///////////

      water_in.f + riserwater_in.f*riserwater_in.z_water - dwater_out.f = der(rho_w*v_l);

      riserwater_in.f*riserwater_in.z_steam - steam_out.f = der(rho_v*v_v);

      V = v_l + v_v;

     // dwater_out.f = 6*riserwater_in.f*riserwater_in.z_steam;

     dwater_out.f = f;

     //steam_out.f=F;

      /////////ENTHALPY BALANCES//////////

      der(rho_w*v_l*cp_water*dwater_out.T) = water_in.f*cp_water*water_in.T + riserwater_in.f*riserwater_in.z_water*riserwater_in.T*cp_water -
      cp_water*dwater_out.T*dwater_out.f;

      der(rho_v*v_v*cp_steam*steam_out.T) = riserwater_in.f*riserwater_in.z_steam*cp_steam*riserwater_in.T - steam_out.f*steam_out.T*cp_steam;

    initial equation

      v_l=1;

     steam_out.T=riserwater_in.T;
    //steam_out.T=100+273;
      dwater_out.T= water_in.T;
     //dwater_out.T= 87+273;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
              Rectangle(
              extent={{-42,70},{52,-6}},
              lineColor={28,108,200},
              fillColor={0,140,72},
              fillPattern=FillPattern.Solid)}), Diagram(coordinateSystem(
              preserveAspectRatio=false), graphics={Rectangle(
              extent={{-42,70},{52,-6}},
              lineColor={28,108,200},
              fillColor={0,140,72},
              fillPattern=FillPattern.Solid)}));
    end steamdrum;
  end equipments;

  package flowsheets
    block economizer_flowsheet
      equipments.economizer economizer
        annotation (Placement(transformation(extent={{-74,-106},{76,92}})));
      streamsources.watersource watersource(f=6.33)
        annotation (Placement(transformation(extent={{96,-6},{76,14}})));
      boiler.streamsources.dryairsource dryairsource
        annotation (Placement(transformation(extent={{-92,30},{-72,50}})));
    equation
      connect(economizer.dryair_in, dryairsource.dryair) annotation (Line(
            points={{-44,38.54},{-62,38.54},{-62,40},{-81.9,40}}, color={28,108,
              200}));
      connect(economizer.water_in, watersource.water) annotation (Line(points={
              {41.5,2.9},{61.75,2.9},{61.75,4},{85.6,4}}, color={28,108,200}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end economizer_flowsheet;

    block superheater_flowsheet
      equipments.superheater superheater
        annotation (Placement(transformation(extent={{-94,-124},{88,84}})));
      streamsources.watersource watersource(
        Water=0,
        Steam=1,
        T=100 + 273)
        annotation (Placement(transformation(extent={{114,-30},{94,-10}})));
      streamsources.watersource watersource1(f=1)
        annotation (Placement(transformation(extent={{106,-70},{86,-50}})));
      boiler.streamsources.dryairsource dryairsource(T=600 + 273, f=100)
        annotation (Placement(transformation(extent={{-108,34},{-88,54}})));

    equation
      connect(superheater.dryair_in, dryairsource.dryair) annotation (Line(
            points={{-66.7,44.48},{-82,44.48},{-82,44},{-97.9,44}}, color={28,
              108,200}));
      connect(watersource.water, superheater.steam_in)
        annotation (Line(points={{103.6,-20},{60.7,-20}}, color={28,108,200}));
      connect(watersource1.water, superheater.attemperator_water) annotation (
          Line(points={{95.6,-60},{76,-60},{76,-51.2},{37.04,-51.2}}, color={28,
              108,200}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end superheater_flowsheet;

    block riser_flowsheet
      equipments.riser riser
        annotation (Placement(transformation(extent={{-90,-106},{90,94}})));
      streamsources.watersource watersource(T=87 + 273, f=100)
                                            annotation (Placement(transformation(extent={{98,-16},
                {78,4}})));
      boiler.streamsources.dryairsource dryairsource(T=500 + 273)
        annotation (Placement(transformation(extent={{-96,44},{-76,64}})));

    equation
      connect(riser.water_in, watersource.water)
        annotation (Line(points={{45,-6},{87.6,-6}}, color={28,108,200}));
      connect(riser.dryair_in, dryairsource.dryair)
        annotation (Line(points={{-55.8,54},{-85.9,54}}, color={28,108,200}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end riser_flowsheet;

    block steamdrum_flowsheet
      equipments.steamdrum steamdrum
        annotation (Placement(transformation(extent={{-72,-134},{62,74}})));
      streamsources.watersource watersource(T=87 + 273, f=6.33)
        annotation (Placement(transformation(extent={{92,24},{72,44}})));
      equipments.riser riser
        annotation (Placement(transformation(extent={{12,-40},{-60,-88}})));
      boiler.streamsources.dryairsource dryairsource(T=500 + 273)
        annotation (Placement(transformation(extent={{50,-90},{30,-70}})));
    equation
      connect(steamdrum.water_in, watersource.water) annotation (Line(points={{
              36.54,32.4},{56.27,32.4},{56.27,34},{81.6,34}}, color={28,108,200}));
      connect(riser.water_out, steamdrum.riserwater_in) annotation (Line(points=
             {{-2.4,-64},{32,-64},{32,-62},{60,-62},{60,-25.84},{36.54,-25.84}},
            color={28,108,200}));
      connect(riser.dryair_in, dryairsource.dryair) annotation (Line(points={{
              -1.68,-78.4},{15.16,-78.4},{15.16,-80},{39.9,-80}}, color={28,108,
              200}));
      connect(steamdrum.dwater_out, riser.water_in) annotation (Line(points={{
              -41.18,-25.84},{-72,-25.84},{-72,-64},{-42,-64}}, color={28,108,
              200}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end steamdrum_flowsheet;

    block overall
      equipments.superheater superheater
        annotation (Placement(transformation(extent={{-88,-12},{-40,62}})));
      equipments.economizer economizer
        annotation (Placement(transformation(extent={{38,-34},{102,72}})));
      equipments.steamdrum steamdrum annotation (Placement(transformation(
            extent={{-30,49},{30,-49}},
            rotation=90,
            origin={20,-51})));
      equipments.riser riser
        annotation (Placement(transformation(extent={{-42,-10},{44,58}})));
      boiler.streamsources.dryairsource dryairsource(T=600 + 273)
        annotation (Placement(transformation(extent={{-114,42},{-94,62}})));
      streamsources.watersource watersource(f=6.33)
        annotation (Placement(transformation(extent={{118,2},{98,22}})));
      streamsources.watersource watersource1
        annotation (Placement(transformation(extent={{-82,-20},{-62,0}})));
    equation
      connect(watersource1.water, superheater.attemperator_water) annotation (
          Line(points={{-71.6,-10},{-40,-10},{-40,13.9},{-53.44,13.9}}, color={
              28,108,200}));
      connect(watersource.water, economizer.water_in) annotation (Line(points={
              {107.6,12},{98,12},{98,24.3},{87.28,24.3}}, color={28,108,200}));
      connect(dryairsource.dryair, superheater.dryair_in) annotation (Line(
            points={{-103.9,52},{-92,52},{-92,47.94},{-80.8,47.94}}, color={28,
              108,200}));
      connect(superheater.dryair_out, riser.dryair_in) annotation (Line(points=
              {{-47.2,47.2},{-37.6,47.2},{-37.6,44.4},{-25.66,44.4}}, color={28,
              108,200}));
      connect(riser.dryair_out, economizer.dryair_in) annotation (Line(points={
              {22.5,44.4},{35.25,44.4},{35.25,43.38},{50.8,43.38}}, color={28,
              108,200}));
      connect(riser.water_out, steamdrum.riserwater_in) annotation (Line(points=
             {{-24.8,24},{-30,24},{-30,-16},{21.96,-16},{21.96,-32.4}}, color={
              28,108,200}));
      connect(economizer.water_out, steamdrum.water_in) annotation (Line(points=
             {{48.88,24.3},{42,24.3},{42,-26},{49.4,-26},{49.4,-32.4}}, color={
              28,108,200}));
      connect(superheater.steam_in, steamdrum.steam_out) annotation (Line(
            points={{-47.2,25},{-36,25},{-36,-86},{49.4,-86},{49.4,-67.2}},
            color={28,108,200}));
      connect(steamdrum.dwater_out, riser.water_in) annotation (Line(points={{
              21.96,-67.2},{21.96,-84},{38,-84},{38,-10},{32,-10},{32,24},{22.5,
              24}}, color={28,108,200}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false), graphics={
            Text(
              extent={{60,16},{80,10}},
              lineColor={28,108,200},
              textString="economizer"),
            Text(
              extent={{-14,12},{6,8}},
              lineColor={28,108,200},
              textString="riser"),
            Text(
              extent={{54,-48},{72,-52}},
              lineColor={28,108,200},
              textString="steam drum"),
            Text(
              extent={{-78,62},{-58,60}},
              lineColor={28,108,200},
              textString="superheater"),
            Text(
              extent={{-116,76},{-102,72}},
              lineColor={28,108,200},
              textString="dry air source"),
            Text(
              extent={{96,-2},{104,-6}},
              lineColor={28,108,200},
              textString="water source"),
            Text(
              extent={{52,-76},{68,-84}},
              lineColor={28,108,200},
              textString="steam out"),
            Text(
              extent={{-2,-70},{14,-76}},
              lineColor={28,108,200},
              textString="downcomer water out")}));
    end overall;

    block boiler_model
      equipments.economizer economizer
        annotation (Placement(transformation(extent={{-26,12},{-104,-70}})));
      equipments.riser riser
        annotation (Placement(transformation(extent={{22,-2},{-30,-66}})));
      equipments.superheater superheater
        annotation (Placement(transformation(extent={{96,6},{32,-66}})));
      streamsources.watersource watersource(f=1)
        annotation (Placement(transformation(extent={{-136,-42},{-116,-22}})));
      boiler.streamsources.dryairsource dryairsource(                   f=10, T=200 +
            273)
        annotation (Placement(transformation(extent={{130,-66},{110,-46}})));
      streamsources.watersource watersource1(f=5)
        annotation (Placement(transformation(extent={{38,-4},{48,6}})));
      equipments.steamdrum steamdrum
        annotation (Placement(transformation(extent={{-44,-46},{32,56}})));
    equation
      connect(dryairsource.dryair, superheater.dryair_in) annotation (Line(
            points={{119.9,-56},{104,-56},{104,-52.32},{86.4,-52.32}}, color={
              28,108,200}));
      connect(watersource1.water, superheater.attemperator_water) annotation (
          Line(points={{43.2,1},{50,1},{50,-19.2},{49.92,-19.2}}, color={28,108,
              200}));
      connect(watersource.water, economizer.water_in) annotation (Line(points={
              {-125.6,-32},{-106,-32},{-106,-33.1},{-86.06,-33.1}}, color={28,
              108,200}));
      connect(superheater.dryair_out, riser.dryair_in) annotation (Line(points=
              {{41.6,-51.6},{26.8,-51.6},{26.8,-53.2},{12.12,-53.2}}, color={28,
              108,200}));
      connect(riser.dryair_out, economizer.dryair_in) annotation (Line(points={
              {-17,-53.2},{-28.5,-53.2},{-28.5,-47.86},{-41.6,-47.86}}, color={
              28,108,200}));
      connect(riser.water_out, steamdrum.riserwater_in) annotation (Line(points=
             {{11.6,-34},{24,-34},{24,7.04},{17.56,7.04}}, color={28,108,200}));
      connect(steamdrum.dwater_out, riser.water_in) annotation (Line(points={{
              -26.52,7.04},{-32,7.04},{-32,-34},{-17,-34}}, color={28,108,200}));
      connect(economizer.water_out, steamdrum.water_in) annotation (Line(points=
             {{-39.26,-33.1},{-39.26,54},{24,54},{24,35.6},{17.56,35.6}}, color=
             {28,108,200}));
      connect(steamdrum.steam_out, superheater.steam_in) annotation (Line(
            points={{-26.52,35.6},{-26.52,66},{30,66},{30,-30},{41.6,-30}},
            color={28,108,200}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false), graphics={
            Text(
              extent={{-20,0},{12,-12}},
              lineColor={28,108,200},
              textString="Steam Drum"),
            Text(
              extent={{-78,-56},{-44,-66}},
              lineColor={28,108,200},
              textString="Economizer"),
            Text(
              extent={{-26,-60},{20,-70}},
              lineColor={28,108,200},
              textString="Riser"),
            Text(
              extent={{44,-60},{88,-70}},
              lineColor={28,108,200},
              textString="Superheater"),
            Text(
              extent={{-144,-48},{-106,-50}},
              lineColor={28,108,200},
              textString="Water_source"),
            Text(
              extent={{102,-68},{136,-76}},
              lineColor={28,108,200},
              textString="Dryair_source"),
            Text(
              extent={{54,-6},{86,12}},
              lineColor={28,108,200},
              textString="Attemeperator")}));
    end boiler_model;
  end flowsheets;
end boiler;
