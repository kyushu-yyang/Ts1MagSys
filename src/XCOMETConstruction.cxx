#include <iostream>
#include <cmath>
#include <ctime>
#include <TString.h>
#include "XMatCopper.hpp"
#include "XMatAluminium.hpp"
#include "XMatNbTi.hpp"

#include "XCoilConductor.hpp"
#include "XCoilStrip.hpp"
#include "XCoilShell.hpp"

#include "XRootOutput.hpp"
#include "XQuenchExcept.hpp"
#include "XQuenchLogger.hpp"
#include "XQuenchOutput.hpp"

#include "XCOMETConstruction.h"


XCOMETConstruction :: XCOMETConstruction()
    : fFld(NULL),
      fTs1a(NULL),
      fTs1b(NULL),
      XQuenchTransient()
{
  if (!fFld)  fFld = new XFieldHandle();
  ConstructField(fFld);
}


XCOMETConstruction :: ~XCOMETConstruction()
{
  if (fFld)   delete fFld;
  if (fTs1a)  delete fTs1a;
  if (fTs1b)  delete fTs1b;
}


void XCOMETConstruction ::ConstructTs1a()
{
  const double r = (250.+16./2)*mm;
  const std::string name = "TS1a";
  
  XCoilHandle* coil = new XCoilHandle();
  coil->SetName(name);
  coil->SetCoilSize(0., 2.*M_PI*r, 0.);
  coil->SetMesh(40, 4, 3);
  coil->SetCoilLayers(1);
  coil->SetCoilTurns(40);
  coil->SetMaterialRatio(7.3, 1., 0.9);

  // set coil structure
  coil->AddLayer(1, kStrip, GetStrip(), kAdiabatic, 0.*cm);
  coil->AddLayer(2, kConductor, GetConductor(), kAdiabatic, 0.*mm);
  coil->AddLayer(3, kShell, GetShell(), kAdiabatic, 0.*mm);

  // set processor
  XProcessManager* pro = new XProcessManager();
  pro->SetCoilHandler(coil);
  pro->Initialize();
  pro->SetNbTiIc(14.10e+3);
  
  // get magnetic field map
  fFld->SetTarget(name);
  fFld->Run();
  pro->SetFieldHandler(fFld);
  //pro->SetUniformField(3.);
  //pro->SetUniformField(2.6);
  
  // uniform RRR and magnetic field
  //pro->SetUniformField(2.0);
  pro->SetUniformRRR(kConductor, 400.);
  pro->SetUniformRRR(kStrip, 2000.);
  //pro->SetUniformRRR(kShell, 10.);

  // fill the material info into the coil solver
  if (!fTs1a) fTs1a = new XThermalSolver();
  fTs1a->SetProcessHandle(pro);

  // write geometry
  XQuenchOutput* geo = new XQuenchOutput( "geoTs1a.geo", iOfstream );
  geo->WriteGeometry(pro);
  geo->Close();

  std::cout << " -- TS1a construction finished." << std::endl;
}


void XCOMETConstruction ::ConstructTs1b()
{
  const double r = (250.+48./2)*mm;
  const std::string name = "TS1b";
  
  XCoilHandle* coil = new XCoilHandle();
  coil->SetName(name);
  coil->SetCoilSize(0., 2.*M_PI*r, 0.);
  coil->SetMesh(48, 4, 5);
  coil->SetCoilLayers(3);
  coil->SetCoilTurns(48);
  coil->SetMaterialRatio(7.3, 1., 0.9);

  // set coil structure
  coil->AddLayer(1, kStrip, GetStrip(), kAdiabatic, 0.*cm);
  coil->AddLayer(2, kConductor, GetConductor(), kAdiabatic, 0.*mm);
  coil->AddLayer(3, kConductor, GetConductor(), kAdiabatic, 0.*mm);
  coil->AddLayer(4, kConductor, GetConductor(), kAdiabatic, 0.*mm);
  coil->AddLayer(5, kShell, GetShell(), kAdiabatic, 0.*mm);

  // set processor
  XProcessManager* pro = new XProcessManager();
  pro->SetCoilHandler(coil);
  pro->Initialize();
  pro->SetNbTiIc(14.10e+3);
  
  // get magnetic field map
  fFld->SetTarget(name);
  fFld->Run();
  pro->SetFieldHandler(fFld);
  //pro->SetUniformField(2.6);
  
  // uniform RRR and magnetic field
  pro->SetUniformRRR(kConductor, 400.);
  pro->SetUniformRRR(kStrip, 2000.);
  //pro->SetUniformRRR(kShell, 10.);

  // fill the material info into the coil solver
  if (!fTs1b) fTs1b = new XThermalSolver();
  fTs1b->SetProcessHandle(pro);

  // write geometry
  XQuenchOutput* geo = new XQuenchOutput( "geoTs1b.geo", iOfstream );
  geo->WriteGeometry(pro);
  geo->Close();

  std::cout << " -- TS1b construction finished." << std::endl;
}

void XCOMETConstruction :: ConstructField(XFieldHandle* fld)
{
  fld->SetCurrent(2700.*Amp);
  fld->AddCoil( "CS0", 857.88*mm, 1038.12*mm, 672.*mm, 823.65*mm );
  fld->SetMesh( "CS0", 35, 9 );

  fld->AddCoil( "CS1", -595.25*mm, 795.25*mm, 672.*mm, 823.65*mm );
  fld->SetMesh( "CS1", 90, 9 );

  fld->AddCoil( "MS1", -2121.375*mm, -653.625*mm, 672.*mm, 756.25*mm );
  fld->SetMesh( "MS1", 57, 5 );

  fld->AddCoil( "MS2", -2910.5*mm, -2189.5*mm, 672.*mm, 789.95*mm );
  fld->SetMesh( "MS2", 70, 7 );

  fld->AddCoil( "TS1a", -2930.*mm, -2730.*mm, 250.*mm, 266.*mm );
  fld->SetMesh( "TS1a", 40, 1 );

  fld->AddCoil( "TS1b", -3200.*mm, -2960.*mm, 250.*mm, 298.*mm );
  fld->SetMesh( "TS1b", 48, 3 );

  fld->AddCoil( "TS1c", -3450.*mm, -3250.*mm, 250.*mm, 314.*mm );
  fld->SetMesh( "TS1c", 40, 4 );
}


XCoilBase* XCOMETConstruction :: GetConductor()
{
  XCoilConductor* cdt = new XCoilConductor();
  cdt->SetDimension( 4.73*mm, 15.*mm );
  cdt->SetInsSize( 0.15*mm, 0.15*mm );

  return dynamic_cast<XCoilBase*>(cdt);
}


XCoilBase* XCOMETConstruction :: GetFirstStrip()
{
  XCoilStrip* strip = new XCoilStrip();
  strip->SetDimension( 4.73*mm+0.15*2*mm, 1.*mm );
  strip->SetInsSize( 0., 0.25/2.*mm );

  return dynamic_cast<XCoilBase*>(strip);
}


XCoilBase* XCOMETConstruction :: GetStrip()
{
  XCoilStrip* strip = new XCoilStrip();
  strip->SetDimension( 4.73*mm+0.15*2*mm, 1.*mm );
  strip->SetInsSize( 0., 0.25*mm );

  return dynamic_cast<XCoilBase*>(strip);
}


XCoilBase* XCOMETConstruction :: GetShell()
{
  XCoilShell * shell = new XCoilShell();
  shell->SetDimension( 4.73*mm+0.15*2*mm, 12.*mm );
  shell->SetInsSize( 0., 5.*mm );

  return dynamic_cast<XCoilBase*>(shell);
}


void XCOMETConstruction :: SetQuenchHeating(XThermalSolver* solve)
{
  //const double A = 30e-3 * 4.25e-3 * (1173./400.);
  //const double V = (1229.+45./2) * 2 * M_PI * 1e-3 * A;
  //std::cout << 6.4*4./0.0015 << "    " << 6.4 / V << std::endl;

  double heat = solve->GetProcess()->GetMaterialEntry(solve->GetProcess()->Id(fHotZ,fHotPhi,fHotR))->GetHeat();
  const double mshz = solve->GetProcess()->GetMesh(iZ);
  //if (solve->GetProcess()->GetMaterialEntry(solve->GetProcess()->Id(fHotZ,fHotPhi,fHotR))->GetStatus()!=kNormal)
  solve->GetProcess()->GetMaterialEntry(solve->GetProcess()->Id(fHotZ,fHotPhi,fHotR))->SetHeat( 2.2/(4.7e-3*15e-3*0.26*M_PI*2./4.) + heat );
}


double XCOMETConstruction :: GetCoilResistance(XThermalSolver* solve)
{
  double res = 0.;
  int idz = 0;
  int idp = 0;
  int idr = 0;

  const int mshz = solve->GetProcess()->GetMesh(iZ);
  const int mshp = solve->GetProcess()->GetMesh(iPhi);
  const int mshr = solve->GetProcess()->GetMesh(iR);

  for (unsigned int i=0; i<solve->GetProcess()->GetMaterialEntries(); i++) {
    idz = solve->GetProcess()->GetDimensionEntry(i)->GetId(iZ);
    idp = solve->GetProcess()->GetDimensionEntry(i)->GetId(iPhi);
    idr = solve->GetProcess()->GetDimensionEntry(i)->GetId(iR);

    if (solve->GetProcess()->GetDimensionEntry(i)->GetGeometry()==kConductor &&
        idz>0 && idz<mshz+1 &&
        idp>0 && idp<mshp+1 &&
        idr>0 && idr<mshr+1 )
      res += solve->GetProcess()->GetMaterialEntry(i)->GetResistance();
  }

  return res;
}


void XCOMETConstruction :: UpdateQuench2(XThermalSolver* solve, const double time)
{
  XMatCopper    cu;
  XMatAluminium al;
  XMatNbTi     *sc = new XMatNbTi;

  sc->SetIcAt5Tesla(14.10e+3);

  double T     = 0.;
  double Tcs   = 0.;
  double RRR   = 0.;
  double B     = 0.;
  double R_Al  = 0.;
  double R_Cu  = 0.;
  double Ic    = 0.;
  double Q     = 0.;
  double V     = 0.;
  double R     = 0.;
  double I_Al  = 0.;
  double I_Cu  = 0.;
  double I_sc  = 0.;

  int idz = 0;
  int idp = 0;
  int idr = 0;

  const double Ec = 0.1 * 1e-6;
  const int mshz = solve->GetProcess()->GetMesh(iZ);
  const int mshp = solve->GetProcess()->GetMesh(iPhi);
  const int mshr = solve->GetProcess()->GetMesh(iR);

  const double factor = solve->GetProcess()->GetCoilHandler()->GetApproachZ();

  const double l_Phi  = solve->GetProcess()->GetCoilHandler()->GetCoilSize(iPhi) / solve->GetProcess()->GetMesh(iPhi);
  const double r_Cu   = solve->GetProcess()->GetCoilHandler()->GetMaterialRatio(iCopper);
  const double r_Al   = solve->GetProcess()->GetCoilHandler()->GetMaterialRatio(iAluminium);
  const double r_nbti = solve->GetProcess()->GetCoilHandler()->GetMaterialRatio(iNbTi);

  const double A_cdt  = solve->GetProcess()->GetCoilHandler()->GetCoilType(2)->GetArea();
  const double A_Cu   = A_cdt * r_Cu;
  const double A_Al   = A_cdt * r_Al;

  double Volume = solve->GetProcess()->GetCoilHandler()->GetCoilType(2)->GetTotalArea() * l_Phi;

  for (unsigned int i=0; i<solve->GetProcess()->GetMaterialEntries(); i++) {
    idz = solve->GetProcess()->GetDimensionEntry(i)->GetId(iZ);
    idp = solve->GetProcess()->GetDimensionEntry(i)->GetId(iPhi);
    idr = solve->GetProcess()->GetDimensionEntry(i)->GetId(iR);

    if ( solve->GetProcess()->GetDimensionEntry(i)->GetGeometry()==kConductor &&
         idz>0 && idz<mshz+1 &&
         idp>0 && idp<mshp+1 &&
         idr>0 && idr<mshr+1 ) {
      T   = solve->GetProcess()->GetMaterialEntry(i)->GetTemperature();
      RRR = solve->GetProcess()->GetMaterialEntry(i)->GetRRR();
      B   = solve->GetProcess()->GetMaterialEntry(i)->GetField();

      // update Temperature, RRR and Field
      cu.SetMaterialProperty(T, 50., B);
      al.SetMaterialProperty(T, RRR, B);
      sc->SetField(B);
      sc->SetTemperature(T);

      Ic  = sc->GetCriticalI();
      //Tcs = sc->GetSharingT( fCurr, 4.5 );

      // calculate resistance
      R_Al  = al.hust_eq_resist(T,RRR,B) * l_Phi / A_Al;
      R_Cu  = cu.GetResistivity() * l_Phi / A_Cu;
      R    = R_Al * R_Cu / (R_Al + R_Cu);

      // check quench states
      if ( fCurr<Ic && Ic>0. ) {
        I_Al = 0.;
        I_Cu = 0.;
        R    = 0.;
        V    = l_Phi * Ec * pow(fCurr/Ic,1.7); 
        if (V/l_Phi<0.01) V = 0.;
        Q    = V * fCurr / Volume;
        solve->GetProcess()->GetMaterialEntry(i)->SetStatus(kSuperconduct);
      }
      else if ( fCurr>=Ic && Ic>0. ) {
        I_sc = Ic;
        I_Al = (fCurr - I_sc) * R_Cu / (R_Al + R_Cu);
        I_Cu = fCurr - I_Al - I_sc;
        //if ( I_Cu<0. ) std::cout << "Error: I_Cu<0.. --> I_Cu: " << I_Cu << std::endl;
        V    = I_Al * R_Al;
        Q    = (R_Al*pow(I_Al,2) + R_Cu*I_Cu*(fCurr-I_Al)) / Volume;
        solve->GetProcess()->GetMaterialEntry(i)->SetStatus(kNormal);
        if ( solve->GetProcess()->GetMaterialEntry(i)->GetQuenchTime()<0. )
          solve->GetProcess()->GetMaterialEntry(i)->SetQuenchTime(time);
      }
      else if ( Ic<=0. ) {
        I_sc = 0.;
        I_Al = (fCurr - I_sc) * R_Cu / (R_Al + R_Cu);
        I_Cu = fCurr - I_Al - I_sc;
        //if ( I_Cu<0. ) std::cout << "Error: I_Cu<0.. --> I_Cu: " << I_Cu << std::endl;
        V    = I_Al * R_Al;
        Q    = (R_Al*pow(I_Al,2) + R_Cu*pow(I_Cu,2)) / Volume;
        solve->GetProcess()->GetMaterialEntry(i)->SetStatus(kNormal);
      }
      else
        std::cout << "Error: Ic is out of range, Ic:" << Ic << " T:" << T << "B: " << B << std::endl;

      solve->GetProcess()->GetMaterialEntry(i)->SetResistance( factor*R );
      solve->GetProcess()->GetMaterialEntry(i)->SetVoltage( factor*V );
      solve->GetProcess()->GetMaterialEntry(i)->SetHeat( Q );

    }
  }
}


void XCOMETConstruction :: UpdateQuench(XThermalSolver* solve, const double time)
{
  XMatCopper    cu;
  XMatAluminium al;
  XMatNbTi     *sc = new XMatNbTi;

  sc->SetIcAt5Tesla(14.10e+3);

  double T     = 0.;
  double RRR   = 0.;
  double B     = 0.;
  double R_Al  = 0.;
  double R_Cu  = 0.;
  double R_avg = 0.;
  double Tcs   = 0.;
  double Tc    = 0.;
  double Ic    = 0.;
  double Rcs   = 0.;

  int idz = 0;
  int idp = 0;
  int idr = 0;

  const int mshz = solve->GetProcess()->GetMesh(iZ);
  const int mshp = solve->GetProcess()->GetMesh(iPhi);
  const int mshr = solve->GetProcess()->GetMesh(iR);

  const double factor = solve->GetProcess()->GetCoilHandler()->GetApproachZ();

  const double l_Phi = solve->GetProcess()->GetCoilHandler()->GetCoilSize(iPhi) / solve->GetProcess()->GetMesh(iPhi);
  const double r_Cu  = solve->GetProcess()->GetCoilHandler()->GetMaterialRatio(iCopper);
  const double r_Al  = solve->GetProcess()->GetCoilHandler()->GetMaterialRatio(iAluminium);
  //const double A_cdt = solve->GetProcess()->GetCoilHandler()->GetCoilParts(kConductor)->GetArea();
  const double A_cdt = solve->GetProcess()->GetCoilHandler()->GetCoilType(2)->GetArea();

  const double A_Cu  = A_cdt * r_Cu;
  const double A_Al  = A_cdt * r_Al;

  //double Volume = solve->GetProcess()->GetCoilHandler()->GetCoilParts(kConductor)->GetTotalArea() * l_Phi;
  double Volume = solve->GetProcess()->GetCoilHandler()->GetCoilType(2)->GetTotalArea() * l_Phi;

  for (unsigned int i=0; i<solve->GetProcess()->GetMaterialEntries(); i++) {
    idz = solve->GetProcess()->GetDimensionEntry(i)->GetId(iZ);
    idp = solve->GetProcess()->GetDimensionEntry(i)->GetId(iPhi);
    idr = solve->GetProcess()->GetDimensionEntry(i)->GetId(iR);

    if ( solve->GetProcess()->GetDimensionEntry(i)->GetGeometry()==kConductor &&
         idz>0 && idz<mshz+1 &&
         idp>0 && idp<mshp+1 &&
         idr>0 && idr<mshr+1 ) {
      T   = solve->GetProcess()->GetMaterialEntry(i)->GetTemperature();
      RRR = solve->GetProcess()->GetMaterialEntry(i)->GetRRR();
      B   = solve->GetProcess()->GetMaterialEntry(i)->GetField();

      // update Temperature, RRR and Field
      cu.SetMaterialProperty(T, 50., B);
      al.SetMaterialProperty(T, RRR, B);
      sc->SetField(B);
      sc->SetTemperature(T);

      // calculate Tcs and Tc
      Tc  = sc->GetCriticalT();
      Ic  = sc->GetCriticalI();
      if (Ic==0.)
        Tcs = T;
      else
        Tcs = T + (Tc - T)*(1 - fCurr/sc->GetCriticalI());

      // calculate resistance
      R_Al  = al.hust_eq_resist(T,RRR,B) * l_Phi / A_Al;
      R_Cu  = cu.GetResistivity() * l_Phi / A_Cu;
      R_avg = pow(1./R_Al + 1./R_Cu, -1);
      
      // setup quench status
      if ( T<Tcs ) {
        solve->GetProcess()->GetMaterialEntry(i)->SetStatus(kSuperconduct);
        solve->GetProcess()->GetMaterialEntry(i)->SetResistance(0.);
        solve->GetProcess()->GetMaterialEntry(i)->SetVoltage(0.);
        solve->GetProcess()->GetMaterialEntry(i)->SetHeat(0.);
        solve->GetProcess()->GetMaterialEntry(i)->SetTcs(Tcs);
      }
      else if ( T>=Tcs && T<Tc ) {
        Rcs = R_avg * (T-Tcs) / (Tc-Tcs);
        solve->GetProcess()->GetMaterialEntry(i)->SetStatus(kTransition);
        solve->GetProcess()->GetMaterialEntry(i)->SetResistance( Rcs*factor );
        solve->GetProcess()->GetMaterialEntry(i)->SetVoltage( factor*Rcs*fCurr );
        solve->GetProcess()->GetMaterialEntry(i)->SetHeat( pow(fCurr,2)*Rcs/Volume );
        solve->GetProcess()->GetMaterialEntry(i)->SetTcs(Tcs);
        if ( solve->GetProcess()->GetMaterialEntry(i)->GetQuenchTime()<0. )
          solve->GetProcess()->GetMaterialEntry(i)->SetQuenchTime(time);
      }
      else if ( T>=Tc ) {
        solve->GetProcess()->GetMaterialEntry(i)->SetStatus(kNormal);
        solve->GetProcess()->GetMaterialEntry(i)->SetResistance(R_avg*factor);
        solve->GetProcess()->GetMaterialEntry(i)->SetVoltage( R_avg*fCurr*factor );
        solve->GetProcess()->GetMaterialEntry(i)->SetHeat( pow(fCurr,2)*R_avg/Volume );
        solve->GetProcess()->GetMaterialEntry(i)->SetTcs(Tcs);
      }

    }
    else {
      solve->GetProcess()->GetMaterialEntry(i)->SetStatus(kSuperconduct);
      solve->GetProcess()->GetMaterialEntry(i)->SetResistance(0.);
      solve->GetProcess()->GetMaterialEntry(i)->SetVoltage(0.);
      solve->GetProcess()->GetMaterialEntry(i)->SetHeat(0.);
    }
  }
}


void XCOMETConstruction :: Begin()
{
  time_t now = time(0);
  tm* local = localtime(&now);
  std::cout << " -- time: " << asctime(local) << std::endl;

  std::cout << "//////////////////////////////////////////////////////////////////////////////" << std::endl;
  std::cout << "//     Quench Simulation Started " << std::endl;
  std::cout << "//////////////////////////////////////////////////////////////////////////////" << std::endl;

  std::cout << " Noting: this is a dedicated study of superconducting magnets quench on COMET" << std::endl;
  std::cout << "         Pion capture solenoid system. The magnet system includes magnets of " << std::endl;
  std::cout << "         CS, MS, TS1a-f. Please see the log file to check the details of par-" << std::endl;
  std::cout << "         ameters of setup." << std::endl;
  std::cout << " Author: Ye Yang (Kyushu University / KEK)" << std::endl;
  std::cout << "         07, July, 2017" << std::endl;
  std::cout << "" << std::endl;
}


void XCOMETConstruction :: Run()
{
  const double Tcool = 4.5*K;

  double dt = fdt;
  fTs1a->SetTimeInterval(dt);
  fTs1b->SetTimeInterval(dt);

  // multiply factor on thermal conductvity of insulation
  fTs1a->GetProcess()->SetInsFactor(2.);
  fTs1b->GetProcess()->SetInsFactor(2.);

  // initial parameters
  double time     = fTime0;
  int    cnt      = 0;
  int    ocnt     = 0;
  double CoilRes  = 0.;
  double qchtime  = fTimef;
  bool   quenched = false;
  bool   preqch   = false;
  double Q        = 0.;
  double dtmin[2] = {0., 0.};
  int    nqch     = 0.;
  double totqch   = GetTotalConductor(fTs1a) + GetTotalConductor(fTs1b);

  while (time<fTimef) {
    
    // 1. update material thermal conductivity, capacity
    fTs1a->GetProcess()->SetMaterial();
    fTs1b->GetProcess()->SetMaterial();

    // 2. update material resistivity, voltage
    UpdateQuench2(fTs1a, time);
    UpdateQuench2(fTs1b, time);

    // 3. calculate coil resistance
    CoilRes = 0.;
    CoilRes += GetCoilResistance(fTs1a);
    CoilRes += GetCoilResistance(fTs1b);

    preqch = quenched;

    if ( fVth<CoilRes*fCurr )
      quenched = true;
    else
      quenched = false;

    if ( preqch==false && quenched==true )
      qchtime = time + fDetTime;

    // 4. calculate the current decay
    if ( quenched==true && time>qchtime ) {
      fCurr = CalCurrentDecay(fPreI, CoilRes, dt);

    // 5. calculate the field decay
      CalFieldDecay(fTs1a);
      CalFieldDecay(fTs1b);
      fPreI = fCurr;
    }

    // set heat generation before quench
    if ( time<=1. )
      SetQuenchHeating(fTs1a);

    // 6. solve the thermal equation
    fTs1a->Solve(dt);
    fTs1a->SetBoundary();

    fTs1b->Solve(dt);
    fTs1b->SetBoundary();

    // 7. connect the shell to each other
    ConnectShell( fTs1b, fTs1a, 30.*mm, dt);

    // calculate the number of quenched conductor
    nqch = GetQuenchConductor(fTs1a, kNormal) + GetQuenchConductor(fTs1b, kNormal);

    // set print out
    if (dt>0.01) fDisplay=1;

    //if (cnt%fDisplay==0 && (int)(time*10000)%10==0) {
    if (cnt%fDisplay==0 || dt>0.1) {
      Q = fTs1a->GetProcess()->GetMaterialEntry(fTs1a->GetProcess()->Id(fHotZ,fHotPhi,fHotR))->GetHeat()*1e-3;
      std::cout << std::setprecision(5)     << "time: " << time << " [sec], Rqch: " 
                << (double)nqch/totqch*100. << " [%], Rtot: "
                << CoilRes << " [Ohm], Vtot: " << CoilRes*fCurr << " [V], I: "
                << fCurr   << " [A], Q: " << Q << " [kW/m3],";
      fTs1a->Print(fHotZ,fHotPhi,fHotR);
    }

    if (cnt%(fDisplay*50)==0 || dt>0.1) {
      XRootOutput output;
      if (ocnt==0)
        output.SetPath("./output");
      output.SetFilename( Form("./output/qchout%i.root",ocnt) );
      output.SetSubDirectory("TS1a");
      output.SetSubDirectory("TS1b");
      output.SetHeader(cnt, time, fCurr, CoilRes, CoilRes*fCurr);
      output.Fill("TS1a", fTs1a->GetProcess());
      output.Fill("TS1b", fTs1b->GetProcess());
      output.Close();
      ocnt++;
    }

    // find the minimum time step
    dtmin[0] = fTs1a->FindTimeStep();
    dtmin[1] = fTs1b->FindTimeStep();

    dt = dtmin[0];
    for (int i=1; i<2; i++) {
      //dt = dtmin[i]<dt ? dtmin[i] : dt;
      if (dtmin[i]<dt) dt = dtmin[i];
    }
    
    time += dt;
    cnt ++;
  }

}


void XCOMETConstruction :: End()
{
  std::cout << "//////////////////////////////////////////////////////////////////////////////" << std::endl;
  std::cout << "//          Quench Simulation Finished" << std::endl;
  std::cout << "//////////////////////////////////////////////////////////////////////////////" << std::endl;

  time_t now = time(0);
  tm* local = localtime(&now);
  std::cout << " -- time: " << asctime(local) << std::endl;
}

void XCOMETConstruction :: ConnectShell(XThermalSolver* mag1, XThermalSolver* mag2, const double l, const double dt)
{
  /*
   * configuration:
         magnet 1               magnet 2
     ||||||||||||||||*~~~~~~~~*|||||||||||           <- shell
     ================          ===========           <- coil
     ================          ===========
     ================          ===========
   */
  
  const int z1    = mag1->GetProcess()->GetMesh(iZ);
  const int r1    = mag1->GetProcess()->GetMesh(iR);
  const int mshp1 = mag1->GetProcess()->GetMesh(iPhi);

  const int z2    = 1;
  const int r2    = mag2->GetProcess()->GetMesh(iR);
  const int mshp2 = mag2->GetProcess()->GetMesh(iPhi);

  if (mshp1!=mshp2) 
    throw;

  int    id   = 0;
  double k    = 0.;
  double C    = 0.;
  double rho  = 0.;
  double alp  = 0.;
  double T1   = 0.;
  double T2   = 0.;
  double x1   = 0.;
  double x2   = 0.;
  double preT = 0.;
  double T    = 0.;

  for (int j=0; j<mshp1; j++) {
    id  = mag1->GetProcess()->Id( z1+1, j+1, r1 );
    preT= mag1->GetProcess()->GetMaterialEntry(id)->GetTemperature();

    // set 1st coil
    id  = mag1->GetProcess()->Id( z1, j+1, r1 );
    T1  = mag1->GetProcess()->GetMaterialEntry(id)->GetTemperature();
    k   = mag1->GetProcess()->GetMaterialEntry(id)->GetConductivity(iZ);
    C   = mag1->GetProcess()->GetMaterialEntry(id)->GetCapacity();
    rho = mag1->GetProcess()->GetMaterialEntry(id)->GetDensity();
    x1  = mag1->GetProcess()->GetDimensionEntry(id)->GetCellSize(iZ) / 2. + l/2.;
    alp = k / rho / C;

    // set 2nd coil
    id  = mag2->GetProcess()->Id( z2, j+1, r2 );
    T2  = mag2->GetProcess()->GetMaterialEntry(id)->GetTemperature();
    x2  = mag2->GetProcess()->GetDimensionEntry(id)->GetCellSize(iZ) / 2. + l/2.;
    
    T   = preT + alp * ((T1-preT)/x1/l - (preT-T2)/x2/l) * dt;

    std::cout << " phi: " << j+1 << ", T: " << T << ", T1: " << T1 << ", T2: " << T2 << std::endl;
    //          << " tmin: " << 1./(alp/x1/l)/2. << ", dt: " << dt << std::endl; 

    id  = mag1->GetProcess()->Id( z1+1, j+1, r1 );
    mag1->GetProcess()->GetMaterialEntry(id)->SetTemperature(T);

    id  = mag1->GetProcess()->Id( z2-1, j+1, r2 );
    mag2->GetProcess()->GetMaterialEntry(id)->SetTemperature(T);
  }
}

