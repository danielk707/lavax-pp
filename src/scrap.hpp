// namespace Eigen {
// template<> struct NumTraits<units::length::angstrom_t>
//  : NumTraits<double> // permits to get the epsilon, dummy_precision, lowest, highest functions
// {
//   typedef units::length::angstrom_t Real;
//   typedef units::length::angstrom_t NonInteger;
//   typedef units::length::angstrom_t Nested;
//   enum {
//     IsComplex = 0,
//     IsInteger = 0,
//     IsSigned = 1,
//     RequireInitialization = 1,
//     ReadCost = 1,
//     AddCost = 3,
//     MulCost = 3
//   };
// };

//   template<> struct NumTraits<boost::units::quantity<boost::units::si::length>>
//  : NumTraits<double> // permits to get the epsilon, dummy_precision, lowest, highest functions
// {
//   typedef units::length::angstrom_t Real;
//   typedef units::length::angstrom_t NonInteger;
//   typedef units::length::angstrom_t Nested;
//   enum {
//     IsComplex = 0,
//     IsInteger = 0,
//     IsSigned = 1,
//     RequireInitialization = 1,
//     ReadCost = 1,
//     AddCost = 3,
//     MulCost = 3
//   };
// };
// }
