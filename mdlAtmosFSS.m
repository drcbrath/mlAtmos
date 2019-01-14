
// US units, use with atmos constructor to construct US units based atmos object
const double Re_us = 20895669.;        // 6369000 / 0.3048;                             // (m) radius of the earth
const double GMR_us = 0.01874329530;   // 0.034163195*(1.8*0.3048);                     // (degR/ft) combined gravity and gas constant of dry air on earth
const double H0_us = 0.0;              //                                               // (ft) datum, sea level
const double T0_us = 518.67;           // 288.15*1.8;                                   // (R), SL std temp
const double rho0_us = 0.00237689;     // 1.225*(0.068521766*0.3048*0.3048*0.3048);     // (sl/ft^3), SL std density
const double P0_us = 2116.2166;        // 101325 * (0.3048*0.3048 / 4.4482216152605);   // (lbf/ft^2), SL std pressure
const double a0_us = 1116.36680;       // 340.2686 / 0.3048;                            // (ft/s), speed of sound at SL std temperature

   // defined atmosphere profiles <need to revise alternate day profiles! base on ref Mil 3013? references !>
const std::vector<double> StdDayHk_us({ 0.000, 3352.800, 6096.000, 9753.600, 14325.600, 15544.800, 21640.800, 25862.890 });
const std::vector<double> StdDayTk_us({ 518.670, 389.970, 389.970, 411.570, 487.170, 487.170, 386.370, 336.510 });
const std::vector<double> StdDayTgradk_us({ 1701.673, 1279.429, 1279.429, 1350.295, 1598.327, 1598.327, 1267.618, 1104.035 });
