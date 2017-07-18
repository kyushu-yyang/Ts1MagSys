/**
 * @file   XCOMETConstruction.h
 * @author Y.Yang (Kyushu University)
 * @date   16 Jul 2017
 **/

#ifndef XCOMETConstruction_HH

#include "XThermalSolver.hpp"
#include "XFieldHandle.hpp"
#include "XCoilHandle.hpp"
#include "XQuenchTransient.hpp"

using namespace Quench;

class XCOMETConstruction : public XQuenchTransient
{
  public:
    /// constructor
    XCOMETConstruction();

    /// deconstructor
    ~XCOMETConstruction();

    /// @brief calculate the coil total resistance
    double GetCoilResistance(XThermalSolver* solver);

    /// @brief construct central coil
    void ConstructTs1a();
    void ConstructTs1b();

    /// @brief update the resistance of conductor
    /// @detail please use this function after SetMaterial() and the current decay calculation
    void UpdateQuench(XThermalSolver* solve, const double time);
    void UpdateQuench2(XThermalSolver* solve, const double time);

    /// @brief begin of run
    virtual void Begin();

    /// @brief run
    virtual void Run();

    /// @brief end of run
    virtual void End();


  protected:
    /// @brief setup conductor size
    XCoilBase* GetConductor();

    /// @brief setup strip size
    XCoilBase* GetStrip();
    XCoilBase* GetFirstStrip();

    /// @brief setup shell size
    XCoilBase* GetShell();

    /// @brief construct magnetic field
    void ConstructField(XFieldHandle* fld);
    void SetQuenchHeating(XThermalSolver* solve);
    
    /// @brief connection of shell
    void ConnectShell(XThermalSolver* mag1, XThermalSolver* mag2, const double l, const double dt);
    void ConnectMagnet(XThermalSolver* mag1, XThermalSolver* mag2, const double l, const double dt);
    

  private:
    XFieldHandle* fFld;
    XThermalSolver* fTs1a;
    XThermalSolver* fTs1b;
};

#endif
