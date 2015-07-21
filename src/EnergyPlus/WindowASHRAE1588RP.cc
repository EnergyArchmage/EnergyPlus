// C++ Headers
#include <string>
#include <iostream>
#include <fstream>

// ObjexxFCL Headers
#include <ObjexxFCL/gio.hh>

// EnergyPlus Headers
#include <WindowASHRAE1588RP.hh>
#include <ConvectionCoefficients.hh>
#include <DataEnvironment.hh>
#include <DataErrorTracking.hh>
#include <DataGlobals.hh>
#include <DataHeatBalance.hh>
#include <DataHeatBalFanSys.hh>
#include <DataHeatBalSurface.hh>
#include <DataIPShortCuts.hh>
#include <DataSurfaces.hh>
#include <DataSystemVariables.hh>
#include <DisplayRoutines.hh>
#include <DataTimings.hh>
#include <General.hh>
#include <HeatBalanceManager.hh>
#include <HeatBalanceSurfaceManager.hh>
#include <InputProcessor.hh>
#include <WindowManager.hh>
#include <SolarShading.hh>
#include <UtilityRoutines.hh>

namespace EnergyPlus {

namespace WindowASHRAE1588RP {

using namespace DataEnvironment;
using namespace DataErrorTracking;
using namespace DataGlobals;
using namespace DataHeatBalance;
using namespace DataHeatBalFanSys;
using namespace DataHeatBalSurface;
using namespace DataIPShortCuts;
using namespace DataSurfaces;
using namespace DataSystemVariables;
using namespace DataTimings;
using namespace General;

using ConvectionCoefficients::SetExtConvectionCoeff;
using ConvectionCoefficients::CalcISO15099WindowIntConvCoeff;
using InputProcessor::GetNumObjectsFound;
using InputProcessor::GetObjectItem;
using InputProcessor::VerifyName;
using HeatBalanceManager::SetupSimpleWindowGlazingSystem;
using HeatBalanceSurfaceManager::InitSolarHeatGains;
using WindowManager::CalcWindowHeatBalance;
using WindowManager::InitGlassOpticalCalculations;
using SolarShading::CalcInteriorSolarDistribution;
using SolarShading::ISABSF;

std::string CurrentModuleObject; // to assist in getting input

void
CreateASHRAE1588RPConstructions( int & ConstrNum, bool & ErrorsFound )
{

  bool stand_alone_analysis = true;  // TODO preprocess this for special executable

  // First read location of the 1588 database file
  if ( GetNumObjectsFound( "DatabaseFile:WindowASHRAE1588RP" ) > 1 ) {
    std::cout << "Error: There should be only one DatabaseFile:WindowASHRAE1588RP object." << std::endl;
    ShowFatalError( "There are multiple DatabaseFile:WindowASHRAE1588RP objects." );
  }

  if ( GetNumObjectsFound( "DatabaseFile:WindowASHRAE1588RP" ) < 1 ) {
    std::cout << "Error: There is no DatabaseFile:WindowASHRAE1588RP object. There must be a DatabaseFile:WindowASHRAE1588RP if your input file contains any Construction:WindowASHRAE1588RP objects." << std::endl;
    ShowFatalError( "There is no DatabaseFile:WindowASHRAE1588RP object. There must be a DatabaseFile:WindowASHRAE1588RP if your input file contains any Construction:WindowASHRAE1588RP objects." );
  }

  int database_num_alpha;
  int database_num_numeric;
  Array1D_string database_alphas( 1 ); // Alpha values array
  Array1D< Real64 > database_numerics( 0 ); // Numeric values array
  int IOStat; // IO Status when calling get input subroutine

  GetObjectItem( "DatabaseFile:WindowASHRAE1588RP", 1, database_alphas, database_num_alpha, database_numerics, database_num_numeric, IOStat, lNumericFieldBlanks, lAlphaFieldBlanks, cAlphaFieldNames, cNumericFieldNames );

  std::string db_1588_file_path_input = database_alphas( 1 );
  std::string db_1588_file_path;
  bool exists;

  std::string RoutineName = "WindowASHRAE1588RP";
  CheckForActualFileName( db_1588_file_path_input, exists, db_1588_file_path );
  if ( ! exists ) {
    ShowSevereError( "WindowASHRAE1588RP: Could not locate the ASHRAE 1588 window database, expecting it as file name=" + db_1588_file_path_input );
    ShowContinueError( "Certain run environments require a full path to be included with the file name in the input field." );
    ShowContinueError( "Try again with putting full path and file name in the field." );
    ShowFatalError( "Program terminates due to these conditions." );
  }

  Json::Value root = read_1588_database(db_1588_file_path);

  // Get the thickness keys from the database
  std::vector < std::string > thickness_keys = root["Glazings"].getMemberNames();
  std::sort(thickness_keys.begin(), thickness_keys.end(),
    [](const std::string &a, const std::string &b) -> bool {
      return std::stod(a) < std::stod(b);
    });

  // Get list of Coatings from the database
  std::vector < std::string > coating_keys = root["Glazings"][thickness_keys[0]].getMemberNames();

  // Get list of Tints from the database
  std::vector < std::string > tint_keys = root["Glazings"][thickness_keys[0]][coating_keys[0]].getMemberNames();

  // Get list of Fenestration Types from the database
  std::vector < std::string > type_keys = root["Types"].getMemberNames();

  // Get list of Frames from the database
  std::vector < std::string > frame_keys = root["Frames"].getMemberNames();

  // Get list of Spacers from the database
  std::vector < std::string > spacer_keys = root["Spacers"].getMemberNames();

  // Get list of Gases from the database
  std::vector < std::string > gas_keys = root["Gases"].getMemberNames();

  // Set testing conditions
  Real64 u_indoor_temp = root["Tests"]["U-factor"]["Indoor Temperature"].asDouble();
  Real64 u_outdoor_temp = root["Tests"]["U-factor"]["Outdoor Temperature"].asDouble();
  Real64 u_wind_speed = root["Tests"]["U-factor"]["Wind Speed"].asDouble();
  Real64 u_solar = root["Tests"]["U-factor"]["Solar Incidence"].asDouble();

  Real64 s_indoor_temp = root["Tests"]["SHGC"]["Indoor Temperature"].asDouble();
  Real64 s_outdoor_temp = root["Tests"]["SHGC"]["Outdoor Temperature"].asDouble();
  Real64 s_wind_speed = root["Tests"]["SHGC"]["Wind Speed"].asDouble();
  Real64 s_solar = root["Tests"]["SHGC"]["Solar Incidence"].asDouble();

  int ConstructNumAlpha; // Number of construction alpha names being passed
  int ConstructNumNumeric; // dummy variable for properties being passed
  Array1D_string ConstructAlphas( 8 ); // Construction Alpha names defined
  Array1D< Real64 > ConstructNumerics( 8 ); // Temporary array to transfer construction properties
  bool ErrorInName;
  bool IsBlank;
  int Loop;

  int TotWinASHRAE1588Constructs = GetNumObjectsFound( "Construction:WindowASHRAE1588RP" ); // Number of window constructions based on ASHRAE 1588RP

  CurrentModuleObject = "Construction:WindowASHRAE1588RP";
  for ( Loop = 1; Loop <= TotWinASHRAE1588Constructs; ++Loop ) { // Loop through all WindowASHRAE1588RP constructions.

    //Get the object names for each construction from the input processor
    GetObjectItem( CurrentModuleObject, Loop, ConstructAlphas, ConstructNumAlpha, ConstructNumerics, ConstructNumNumeric, IOStat, lNumericFieldBlanks, lAlphaFieldBlanks, cAlphaFieldNames, cNumericFieldNames );

    ErrorInName = false;
    IsBlank = false;
    VerifyName( ConstructAlphas( 1 ), Construct.Name(), ConstrNum, ErrorInName, IsBlank, trim( CurrentModuleObject ) + " Name" );
    if ( IsBlank ) {
      ErrorsFound = true;
      continue;
    }

    ++ConstrNum;

    int TotMaterialsSave = TotMaterials;

    // Save Materials

    Array1D< MaterialProperties > MaterialSave;
    Array1D< Real64 > NominalRSave;
    MaterialSave.allocate( TotMaterials );
    NominalRSave.allocate( TotMaterials );
    MaterialSave = Material;
    NominalRSave = NominalR;
    Material.deallocate();
    NominalR.deallocate();

    int number_of_gaps;
    int number_of_new_materials;

    Array1D< MaterialProperties > new_materials;
    Array1D< Real64 > new_nominal_R;

    // Read spectral data from database

    // Save Spectral Data
    int TotSpectralDataSave = TotSpectralData;
    Array1D< SpectralDataProperties > SpectralDataSave;
    SpectralDataSave.allocate( TotSpectralData );
    SpectralDataSave = SpectralData;
    SpectralData.deallocate();

    Array1D< SpectralDataProperties > new_spectraldata;
    int num_spectral_datasets;

    // Save Constructions -- The list will be deleted so that the only
    // construction is the one currently being set for any borrowed subroutines
    Array1D< ConstructionData > ConstructSave;
    Array1D< Real64 > NominalRforNominalUCalculationSave;
    Array1D< Real64 > NominalUSave;
    int TotConstructsSave = TotConstructs;

    {
      ConstructSave.allocate( TotConstructs );
      ConstructSave = Construct;
      NominalRforNominalUCalculationSave.allocate( TotConstructs );
      NominalUSave.allocate( TotConstructs );
      NominalRforNominalUCalculationSave = NominalRforNominalUCalculation;
      NominalUSave = NominalU;

      Construct.deallocate();
      NominalRforNominalUCalculation.deallocate();
      NominalU.deallocate();

      Construct.allocate(1);
      NominalRforNominalUCalculation.allocate(1);
      NominalU.allocate(1);

      TotConstructs = 1;
    }


    ConstructionData new_construct;

    // Name
    std::string construction_name = ConstructAlphas( 1 );

    // U-factor
    Real64 target_u_factor;
    bool u_factor_set;

    if ( lNumericFieldBlanks( 1 ) )
    {
      u_factor_set = false;
    }
    else
    {
      u_factor_set = true;
      target_u_factor = ConstructNumerics( 1 );
    }

    // SHGC
    Real64 target_shgc;
    bool shgc_set;

    if ( lNumericFieldBlanks( 2 ) )
    {
      shgc_set = false;
    }
    else
    {
      shgc_set = true;
      target_shgc = ConstructNumerics( 2 );
    }
    // set initial values and set locks as appropriate. If IDF is blank,
    // override with defaults from ASHRAE 1588 RP Database file (TODO)

    // Fenestration Type
    std::string fenestration_type;
    if ( lAlphaFieldBlanks( 2 ) )
    {
      fenestration_type = "FIXED";
    }
    else
    {
      fenestration_type = ConstructAlphas( 2 );
      search_database_keys_for_input(construction_name, "Fenestration Type", type_keys, fenestration_type, ErrorsFound);
    }

    // Number of Panes
    int number_of_panes;
    bool number_of_panes_lock;

    if ( lNumericFieldBlanks( 3 ) )
    {
      number_of_panes_lock = false;
    }
    else
    {
      number_of_panes_lock = true;
      number_of_panes = ConstructNumerics( 3 );
    }
    if ( ! number_of_panes_lock )
    {
      if (u_factor_set)
      {
        // For now use Arasteh method
        if (target_u_factor < (1.4 + 1.7)/2.0) // use average for now since we can't interpolate
        {
          if (shgc_set)
          {
            if ( target_shgc < (0.35 + 0.45)/2 )
            {
              number_of_panes = 3;
            }
            else
            {
              number_of_panes = 2; // vacuum glazings
            }
          }
        }
        else if (target_u_factor < (3.4 + 4.5)/2.0)
        {
          number_of_panes = 2;
        }
        else
        {
          number_of_panes = 1;
        }
      }
      else
      {
        number_of_panes = 2;
      }
    }

    // Glazing Thickness
    Real64 glass_thickness;
    bool glass_thickness_lock;
    if ( lNumericFieldBlanks( 4 ) )
    {
      glass_thickness_lock = false;
    }
    else
    {
      glass_thickness_lock = true;
      glass_thickness = ConstructNumerics( 4 );
    }

    if ( ! glass_thickness_lock )
    {
      glass_thickness = 0.003;
    }

    // Glazing Tint
    std::string glazing_tint;
    bool glazing_tint_lock;

    if ( lAlphaFieldBlanks( 3 ) )
    {
      glazing_tint_lock = false;
    }
    else
    {
      glazing_tint_lock = true;
      glazing_tint = ConstructAlphas( 3 );
      search_database_keys_for_input(construction_name, "Glazing Tint Type", tint_keys, glazing_tint, ErrorsFound);
    }

    if (! glazing_tint_lock )
    {
      glazing_tint = "CLEAR";
    }

    // Glazing Coating
    std::string glazing_coating;
    bool glazing_coating_lock;

    if ( lAlphaFieldBlanks( 4 ) )
    {
      glazing_coating_lock = false;
    }
    else
    {
      glazing_coating_lock = true;
      glazing_coating = ConstructAlphas( 4 );
      search_database_keys_for_input(construction_name, "Glazing Coating Type", coating_keys, glazing_coating, ErrorsFound);
    }

    if (! glazing_coating_lock )
    {
      glazing_coating = "NONE";
    }

    // Gas Type
    std::string gas_type;
    bool gas_type_lock;

    if ( lAlphaFieldBlanks( 5 ) )
    {
      gas_type_lock = false;
    }
    else
    {
      gas_type_lock = true;
      gas_type = ConstructAlphas( 5 );
      search_database_keys_for_input(construction_name, "Gas Type", gas_keys, gas_type, ErrorsFound);
    }

    if (! gas_type_lock )
    {
      gas_type = "AIR";
    }

    // Gap Thickness
    Real64 gap_thickness;
    bool gap_thickness_lock;
    if ( lNumericFieldBlanks( 5 ) )
    {
      gap_thickness_lock = false;
    }
    else
    {
      gap_thickness_lock = true;
      gap_thickness = ConstructNumerics( 5 );
    }

    if ( ! gap_thickness_lock )
    {
      gap_thickness = 0.0127;
    }

    // Spacer Material Type
    std::string spacer_type;
    bool spacer_type_lock;

    if ( lAlphaFieldBlanks( 6 ) )
    {
      spacer_type_lock = false;
    }
    else
    {
      spacer_type_lock = true;
      spacer_type = ConstructAlphas( 6 );
      search_database_keys_for_input(construction_name, "Spacer Material Type", spacer_keys, spacer_type, ErrorsFound);
    }

    if (! spacer_type_lock )
    {
      spacer_type = "STEEL";
    }

    // Frame Material
    std::string frame_material;
    bool frame_material_lock;

    if ( lAlphaFieldBlanks( 7 ) )
    {
      frame_material_lock = false;
    }
    else
    {
      frame_material_lock = true;
      frame_material = ConstructAlphas( 7 );
      search_database_keys_for_input(construction_name, "Frame Material Type", frame_keys, frame_material, ErrorsFound);
    }

    if (! frame_material_lock )
    {
      frame_material = "VINYL";
    }

    // Frame Width
    Real64 frame_width;
    bool frame_width_lock;

    if ( lNumericFieldBlanks( 6 ) )
    {
      frame_width_lock = false;
    }
    else
    {
      frame_width_lock = true;
      frame_width = ConstructNumerics( 6 );
    }

    if (! frame_width_lock )
    {
      frame_width = 0.05;
    }

    // Divider Width
    Real64 divider_width;
    bool divider_width_lock;

    if ( lNumericFieldBlanks( 7 ) )
    {
      divider_width_lock = false;
    }
    else
    {
      divider_width_lock = true;
      divider_width = ConstructNumerics( 7 );
    }

    if (! divider_width_lock )
    {
      divider_width = 0.0;
    }

    // Dirt Factor
    Real64 dirt_factor;
    if ( lNumericFieldBlanks( 8 ) )
    {
      dirt_factor = 1.0;
    }
    else
    {
      dirt_factor = ConstructNumerics( 8 );
    }

    std::string ashrae1588_file_name;
    if ( lAlphaFieldBlanks( 8 ) )
    {
      ashrae1588_file_name = "";
    }
    else
    {
      ashrae1588_file_name = ConstructAlphas(8);
    }

   // Check for Errors
   if ( ErrorsFound ) {
     std::cout << "Error found in processing ASHRAE 1588 window construction input." << std::endl;
     ShowFatalError( "Error found in processing ASHRAE 1588 window construction input." );
   }

    new_construct.Name = construction_name;
    new_construct.TypeIsWindow = true;

    // Save Frame and Divider objects
    int TotFrameDividerSave = TotFrameDivider;
    Array1D< FrameDividerProperties > FrameDividerSave;

    FrameDividerSave.allocate( TotFrameDivider );
    FrameDividerSave = FrameDivider;

    FrameDivider.deallocate();

    FrameDivider.allocate(1);

    TotFrameDivider = 1;

    FrameDividerProperties new_frame_divider;

    Real64 frame_solar_absorptivity;
    Real64 frame_visible_absorptivity;


    // internal defaults based on other values
    Real64 frame_conductance;
    Real64 frame_edge_ratio;

    Real64 fenestration_width;
    Real64 fenestration_height;
    Real64 glazing_width;
    Real64 glazing_height;
    Real64 tilt;
    Real64 fenestration_area;
    Real64 glazing_area;
    Real64 frame_area;

    bool has_frame;
    int num_horizontal_dividers;
    int num_vertical_dividers;

    // matching variables
    Real64 u_factor;
    Real64 u_cog;
    Real64 u_eog;
    Real64 shgc;
    Real64 vt;

    // internal defaults to be left alone
    Real64 frame_IR_emissivity = 0.9;

    Real64 max_divider_spacing = 0.3; // NFRC 100-2014 4.2.2 (B)
    Real64 edge_width = 0.06355;

    // Allocate temporary arrays
    create_dummy_variables();

    Surface( 1 ).Name = construction_name + ":Surface";

    ASHRAE1588RP_Flag = true;
    KickOffSimulation = false;

    Real64 u_factor_match_tolerance = 0.05; // Precision of NFRC reporting TODO expose?
    Real64 optical_match_tolerance = 0.01; // Precision of NFRC reporting TODO expose?

    Real64 u_factor_diff;
    Real64 shgc_diff;

    bool target_matched = false;

    // This is where the iterative optimization loop will begin
    while (! target_matched)
    {

      // Get nearest thickness
      std::string thickness_key;
      if (glass_thickness <= std::stod(thickness_keys[0])/1000.0)
        thickness_key = thickness_keys[0];
      else if (glass_thickness >= std::stod(thickness_keys.back())/1000.0)
        thickness_key = thickness_keys.back();
      else {
        for (std::vector< std::string >::const_iterator i = thickness_keys.begin(); i != thickness_keys.end(); ++i) {
          double first = std::stod(*i)/1000.0;
          double second = std::stod(*std::next(i))/1000.0;
          if (glass_thickness >= first && glass_thickness <= second) {
            if (glass_thickness - first < second - glass_thickness)
              thickness_key = *i;
            else
              thickness_key = *std::next(i);
          }
        }
      }

      frame_solar_absorptivity = root["Types"][fenestration_type]["Frame Absorptivity"].asDouble();
      frame_visible_absorptivity = root["Types"][fenestration_type]["Frame Absorptivity"].asDouble();

      Real64 frame_u_factor;
      std::string fenestration_category = root["Types"][fenestration_type]["Category"].asString();
      if ( root["Spacers"][spacer_type]["Metal"].asBool() ) {
        frame_u_factor = root["Frames"][frame_material]["Metal Spacer"][fenestration_category][std::max(number_of_panes,3)-1].asDouble();
      }
      else {
        frame_u_factor = root["Frames"][frame_material]["Non-Metal Spacer"][fenestration_category][std::max(number_of_panes,3)-1].asDouble();
      }

      Real64 assumed_h_o = 30;  // W/m2-K
      Real64 assumed_h_i = 8;  // W/m2-K
      frame_conductance = 1/(1/frame_u_factor - (1/assumed_h_o + 1/assumed_h_i));

      fenestration_width = root["Types"][fenestration_type]["Width"].asDouble();
      fenestration_height = root["Types"][fenestration_type]["Height"].asDouble();
      tilt = root["Types"][fenestration_type]["Tilt"].asDouble()*Pi/180.0;

      if ( frame_width > 0.0 )
      {
        has_frame = true;
      }
      else
      {
        has_frame = false;
      }

      fenestration_area = fenestration_width*fenestration_height;
      glazing_width = fenestration_width - 2.0*frame_width;
      glazing_height = fenestration_height - 2.0*frame_width;
      glazing_area = glazing_width*glazing_height;

      if ( has_frame && divider_width > 0.0)
      {
        num_horizontal_dividers = ceil(glazing_height/max_divider_spacing);
        num_vertical_dividers = ceil(glazing_width/max_divider_spacing);
      }
      else
      {
        num_horizontal_dividers = 0;
        num_vertical_dividers = 0;
      }


      Surface( 1 ).Height = glazing_height;
      Surface( 1 ).Width = glazing_width;
      Surface( 1 ).Area = glazing_area;
      Surface( 1 ).Tilt = tilt*180/Pi;
      Surface( 1 ).CosTilt = cos(tilt);
      Surface( 1 ).SinTilt = sin(tilt);
      Surface( 1 ).ViewFactorSky = 0.5 * ( 1.0 + Surface( 1 ).CosTilt );
      Surface( 1 ).ViewFactorGround = 0.5 * ( 1.0 - Surface( 1 ).CosTilt );
      Surface( 1 ).ViewFactorSkyIR = Surface( 1 ).ViewFactorSky;
      Surface( 1 ).ViewFactorGroundIR = Surface( 1 ).ViewFactorGround;
      AirSkyRadSplit( 1 ) = std::sqrt( 0.5 * ( 1.0 + Surface( 1 ).CosTilt ) );

      number_of_gaps = number_of_panes - 1;
      number_of_new_materials = number_of_panes + number_of_gaps;

      // Construction specific allocations
      AWinSurf.allocate(number_of_panes, 1);
      QRadSWwinAbs.allocate(number_of_panes, 1);
      QRadSWwinAbsLayer.allocate(number_of_panes, 1);

      // Create New Spectral Data objects
      num_spectral_datasets = std::min(number_of_panes,2);
      if ( new_spectraldata.size_ != (unsigned)num_spectral_datasets ) {
        new_spectraldata.allocate(num_spectral_datasets);
        SpectralData.allocate(num_spectral_datasets);
        TotSpectralData = num_spectral_datasets;
      }

      SpectralData( 1 ).Name = construction_name + ":SPECTRALDATA1";
      int sd1_num_wavelengths = root["Glazings"][thickness_key][glazing_coating][glazing_tint]["Wavelengths"].size();
      SpectralData( 1 ).NumOfWavelengths = sd1_num_wavelengths;

      SpectralData( 1 ).WaveLength.allocate( sd1_num_wavelengths ); // Wavelength (microns)
      SpectralData( 1 ).Trans.allocate( sd1_num_wavelengths ); // Transmittance at normal incidence
      SpectralData( 1 ).ReflFront.allocate( sd1_num_wavelengths ); // Front reflectance at normal incidence
      SpectralData( 1 ).ReflBack.allocate( sd1_num_wavelengths ); // Back reflectance at normal incidence

      for ( int i = 1; i <= sd1_num_wavelengths; i++ ) {
        SpectralData( 1 ).WaveLength( i ) = root["Glazings"][thickness_key][glazing_coating][glazing_tint]["Wavelengths"][i-1].asDouble(); // Wavelengths
        SpectralData( 1 ).Trans( i ) = root["Glazings"][thickness_key][glazing_coating][glazing_tint]["T"][i-1].asDouble();
        // Following is needed since angular calculation in subr TransAndReflAtPhi
        // fails for Trans = 0.0
        if ( SpectralData( 1 ).Trans( i ) < 0.001 ) SpectralData( 1 ).Trans( i ) = 0.001;
        SpectralData( 1 ).ReflFront( i ) = root["Glazings"][thickness_key][glazing_coating][glazing_tint]["Rf"][i-1].asDouble();
        SpectralData( 1 ).ReflBack( i ) = root["Glazings"][thickness_key][glazing_coating][glazing_tint]["Rb"][i-1].asDouble();
      }

      // Check integrity of the spectral data
      for ( int LamNum = 1; LamNum <= sd1_num_wavelengths; ++LamNum ) {
        Real64 Lam = SpectralData( 1 ).WaveLength( LamNum );
        Real64 Tau = SpectralData( 1 ).Trans( LamNum );
        Real64 RhoF = SpectralData( 1 ).ReflFront( LamNum );
        Real64 RhoB = SpectralData( 1 ).ReflBack( LamNum );
        if ( LamNum < sd1_num_wavelengths ) {
          if ( SpectralData( 1 ).WaveLength( LamNum + 1 ) <= Lam ) {
            ErrorsFound = true;
            ShowSevereError( RoutineName + CurrentModuleObject + "=\"" + SpectralData( 1 ).Name + "\" invalid set." );
            ShowContinueError( "... Wavelengths not in increasing order. at wavelength#=" + TrimSigDigits( LamNum ) + ", value=[" + TrimSigDigits( Lam, 4 ) + "], next is [" + TrimSigDigits( SpectralData( Loop ).WaveLength( LamNum + 1 ), 4 ) + "]." );
          }
        }

        if ( Tau > 1.01 ) {
          ErrorsFound = true;
          ShowSevereError( RoutineName + CurrentModuleObject + "=\"" + SpectralData( 1 ).Name + "\" invalid value." );
          ShowContinueError( "... A transmittance is > 1.0; at wavelength#=" + TrimSigDigits( LamNum ) + ", value=[" + TrimSigDigits( Tau, 4 ) + "]." );
        }

        if ( RhoF < 0.0 || RhoF > 1.02 || RhoB < 0.0 || RhoB > 1.02 ) {
          ErrorsFound = true;
          ShowSevereError( RoutineName + CurrentModuleObject + "=\"" + SpectralData( 1 ).Name + "\" invalid value." );
          ShowContinueError( "... A reflectance is < 0.0 or > 1.0; at wavelength#=" + TrimSigDigits( LamNum ) + ", RhoF value=[" + TrimSigDigits( RhoF, 4 ) + "]." );
          ShowContinueError( "... A reflectance is < 0.0 or > 1.0; at wavelength#=" + TrimSigDigits( LamNum ) + ", RhoB value=[" + TrimSigDigits( RhoB, 4 ) + "]." );
        }

        if ( ( Tau + RhoF ) > 1.03 || ( Tau + RhoB ) > 1.03 ) {
          ErrorsFound = true;
          ShowSevereError( RoutineName + CurrentModuleObject + "=\"" + SpectralData( 1 ).Name + "\" invalid value." );
          ShowContinueError( "... Transmittance + reflectance) > 1.0 for an entry; at wavelength#=" + TrimSigDigits( LamNum ) + ", value(Tau+RhoF)=[" + TrimSigDigits( ( Tau + RhoF ), 4 ) + "], value(Tau+RhoB)=[" + TrimSigDigits( ( Tau + RhoB ), 4 ) + "]." );
        }

      }

      if ( num_spectral_datasets == 2 ) {
        SpectralData( 2 ).Name = construction_name + ":SPECTRALDATA2";
        int sd2_num_wavelengths = root["Glazings"][thickness_key]["NONE"]["CLEAR"]["Wavelengths"].size();
        SpectralData( 2 ).NumOfWavelengths = sd2_num_wavelengths;

        SpectralData( 2 ).WaveLength.allocate( sd1_num_wavelengths ); // Wavelength (microns)
        SpectralData( 2 ).Trans.allocate( sd1_num_wavelengths ); // Transmittance at normal incidence
        SpectralData( 2 ).ReflFront.allocate( sd1_num_wavelengths ); // Front reflectance at normal incidence
        SpectralData( 2 ).ReflBack.allocate( sd1_num_wavelengths ); // Back reflectance at normal incidence

        for ( int i = 1; i <= sd1_num_wavelengths; i++ ) {
          SpectralData( 2 ).WaveLength( i ) = root["Glazings"][thickness_key]["NONE"]["CLEAR"]["Wavelengths"][i-1].asDouble();
          SpectralData( 2 ).Trans( i ) = root["Glazings"][thickness_key]["NONE"]["CLEAR"]["T"][i-1].asDouble();
          // Following is needed since angular calculation in subr TransAndReflAtPhi
          // fails for Trans = 0.0
          if ( SpectralData( 2 ).Trans( i ) < 0.001 ) SpectralData( 2 ).Trans( i ) = 0.001;
          SpectralData( 2 ).ReflFront( i ) = root["Glazings"][thickness_key]["NONE"]["CLEAR"]["Rf"][i-1].asDouble();
          SpectralData( 2 ).ReflBack( i ) = root["Glazings"][thickness_key]["NONE"]["CLEAR"]["Rb"][i-1].asDouble();
        }

        // Check integrity of the spectral data
        for ( int LamNum = 1; LamNum <= sd2_num_wavelengths; ++LamNum ) {
          Real64 Lam = SpectralData( 2 ).WaveLength( LamNum );
          Real64 Tau = SpectralData( 2 ).Trans( LamNum );
          Real64 RhoF = SpectralData( 2 ).ReflFront( LamNum );
          Real64 RhoB = SpectralData( 2 ).ReflBack( LamNum );
          if ( LamNum < sd2_num_wavelengths ) {
            if ( SpectralData( 2 ).WaveLength( LamNum + 1 ) <= Lam ) {
              ErrorsFound = true;
              ShowSevereError( RoutineName + CurrentModuleObject + "=\"" + SpectralData( 2 ).Name + "\" invalid set." );
              ShowContinueError( "... Wavelengths not in increasing order. at wavelength#=" + TrimSigDigits( LamNum ) + ", value=[" + TrimSigDigits( Lam, 4 ) + "], next is [" + TrimSigDigits( SpectralData( Loop ).WaveLength( LamNum + 1 ), 4 ) + "]." );
            }
          }

          if ( Tau > 1.01 ) {
            ErrorsFound = true;
            ShowSevereError( RoutineName + CurrentModuleObject + "=\"" + SpectralData( 2 ).Name + "\" invalid value." );
            ShowContinueError( "... A transmittance is > 1.0; at wavelength#=" + TrimSigDigits( LamNum ) + ", value=[" + TrimSigDigits( Tau, 4 ) + "]." );
          }

          if ( RhoF < 0.0 || RhoF > 1.02 || RhoB < 0.0 || RhoB > 1.02 ) {
            ErrorsFound = true;
            ShowSevereError( RoutineName + CurrentModuleObject + "=\"" + SpectralData( 2 ).Name + "\" invalid value." );
            ShowContinueError( "... A reflectance is < 0.0 or > 1.0; at wavelength#=" + TrimSigDigits( LamNum ) + ", RhoF value=[" + TrimSigDigits( RhoF, 4 ) + "]." );
            ShowContinueError( "... A reflectance is < 0.0 or > 1.0; at wavelength#=" + TrimSigDigits( LamNum ) + ", RhoB value=[" + TrimSigDigits( RhoB, 4 ) + "]." );
          }

          if ( ( Tau + RhoF ) > 1.03 || ( Tau + RhoB ) > 1.03 ) {
            ErrorsFound = true;
            ShowSevereError( RoutineName + CurrentModuleObject + "=\"" + SpectralData( 2 ).Name + "\" invalid value." );
            ShowContinueError( "... Transmittance + reflectance) > 1.0 for an entry; at wavelength#=" + TrimSigDigits( LamNum ) + ", value(Tau+RhoF)=[" + TrimSigDigits( ( Tau + RhoF ), 4 ) + "], value(Tau+RhoB)=[" + TrimSigDigits( ( Tau + RhoB ), 4 ) + "]." );
          }

        }
      }

      // Create New Material objects
      if ( new_materials.size_ != (unsigned)number_of_new_materials ) {
        new_materials.allocate( number_of_new_materials );
        Material.allocate( number_of_new_materials );
        NominalR.allocate( number_of_new_materials );
        TotMaterials = number_of_new_materials;
      }


      // Define material properties for glazings
      for ( int MaterNum = 1; MaterNum <= number_of_new_materials; MaterNum += 2 )
      {
        std::string coating, tint;
        if (MaterNum == 1) {
          coating = glazing_coating;
          tint = glazing_tint;
          Material( MaterNum ).GlassSpectralDataPtr = 1;
        }
        else {
          coating = "NONE";
          tint = "CLEAR";
          Material( MaterNum ).GlassSpectralDataPtr = 2;
        }
        Material( MaterNum ).Group = WindowGlass;
        Material( MaterNum ).Name = construction_name + ":GLAZING" + std::to_string(MaterNum);
        Material( MaterNum ).Roughness = VerySmooth;
        Material( MaterNum ).ROnly = true;
        Material( MaterNum ).Thickness = glass_thickness;
        Material( MaterNum ).Trans = root["Glazings"][thickness_key][coating][tint]["Tsol"].asDouble();
        Material( MaterNum ).ReflectSolBeamFront = root["Glazings"][thickness_key][coating][tint]["Rfsol"].asDouble();
        Material( MaterNum ).ReflectSolBeamBack = root["Glazings"][thickness_key][coating][tint]["Rbsol"].asDouble();
        Material( MaterNum ).TransVis = root["Glazings"][thickness_key][coating][tint]["Tvis"].asDouble();
        Material( MaterNum ).ReflectVisBeamFront = root["Glazings"][thickness_key][coating][tint]["Rfvis"].asDouble();
        Material( MaterNum ).ReflectVisBeamBack = root["Glazings"][thickness_key][coating][tint]["Rbvis"].asDouble();
        Material( MaterNum ).TransThermal = 0.0;
        Material( MaterNum ).AbsorpThermalFront = root["Glazings"][thickness_key][coating][tint]["ef"].asDouble();
        Material( MaterNum ).AbsorpThermalBack = root["Glazings"][thickness_key][coating][tint]["eb"].asDouble();
        Material( MaterNum ).Conductivity = 1.0;
        Material( MaterNum ).GlassTransDirtFactor = 1.0;  // Hold at unity to find match and then apply to outside layer
        Material( MaterNum ).YoungModulus = 7.2e10;
        Material( MaterNum ).PoissonsRatio = 0.22;
        Material( MaterNum ).AbsorpThermal = Material( MaterNum ).AbsorpThermalBack;
        Material( MaterNum ).SolarDiffusing = false;

        NominalR( MaterNum ) = Material( MaterNum ).Thickness / Material( MaterNum ).Conductivity;
        Material( MaterNum ).Resistance = NominalR( MaterNum );

      }

      int num_gases = root["Gases"][gas_type].size();

      // Define material properties for gaps
      for ( int MaterNum = 2; MaterNum <= number_of_new_materials; MaterNum += 2 )
      {
        Material( MaterNum ).Group = WindowGas;
        Material( MaterNum ).Name = construction_name + ":GAP" + std::to_string(MaterNum);
        Material( MaterNum ).Roughness = MediumRough;
        Material( MaterNum ).ROnly = true;
        Material( MaterNum ).Thickness = gap_thickness;
        Material( MaterNum ).NumberOfGasesInMixture = num_gases;

        for ( int gas = 1; gas <= num_gases; gas++)
        {
          Material( MaterNum ).GasFract( gas ) = root["Gases"][gas_type][gas-1]["Fraction"].asDouble();
          Material( MaterNum ).GasWght( gas ) = root["Gases"][gas_type][gas-1]["Molecular Weight"].asDouble();
          Material( MaterNum ).GasSpecHeatRatio( gas ) = root["Gases"][gas_type][gas-1]["Specific Heat Ratio"].asDouble();
          for ( int ICoeff = 1; ICoeff <= 3; ++ICoeff ) {
            Material( MaterNum ).GasCon( gas, ICoeff ) = root["Gases"][gas_type][gas-1]["Conductivity"][ICoeff-1].asDouble();
            Material( MaterNum ).GasVis( gas, ICoeff ) = root["Gases"][gas_type][gas-1]["Viscosity"][ICoeff-1].asDouble();
            Material( MaterNum ).GasCp( gas, ICoeff ) = root["Gases"][gas_type][gas-1]["Specific Heat"][ICoeff-1].asDouble();
          }
        }

        Real64 DenomRGas = ( Material( MaterNum ).GasCon( 1, 1 ) + Material( MaterNum ).GasCon( 1, 2 ) * 300.0 + Material( MaterNum ).GasCon( 1, 3 ) * 90000.0 );
        NominalR( MaterNum ) = Material( MaterNum ).Thickness / DenomRGas;

      }

      new_construct.TotLayers = number_of_new_materials;

      for ( int Layer = 1; Layer <= number_of_new_materials; ++Layer ) {
        new_construct.LayerPoint( Layer ) = Layer;
      }

      Construct( 1 ) = new_construct;

      NominalRforNominalUCalculation( 1 ) = 0.0;
      for ( int Layer = 1; Layer <= Construct( 1 ).TotLayers; ++Layer ) {
        NominalRforNominalUCalculation( 1 ) += NominalR( Construct( 1 ).LayerPoint( Layer ) );
      }

      CheckAndSetConstructionProperties( 1, ErrorsFound );

      Surface( 1 ).Construction = 1; // This is the only construction available to the dummy surface. The actual surface will reference the real construction.
      Surface( 1 ).FrameDivider = 0; // Set temporarily until after Center-of-Glass U-factor is calculated

      // Set up U-factor conditions TODO read these from database
      Real64 in_air_temp = u_indoor_temp;
      Real64 out_air_temp = u_outdoor_temp;
      Real64 wind_speed = u_wind_speed;
      Real64 solar_incident = u_solar;

      // Calculate Center-of-Glass U-factor (without Frame)
      calc_window_performance(in_air_temp, out_air_temp, wind_speed, solar_incident);

      u_cog = -WinHeatGain(1)/(Surface( 1 ).Area*(in_air_temp - out_air_temp));

      if ( number_of_panes == 1)
        u_eog = u_cog;
      else {
        Real64 eog_a = root["Spacers"][spacer_type]["Coefficients"][std::max(number_of_panes,3)-2][0].asDouble();
        Real64 eog_b = root["Spacers"][spacer_type]["Coefficients"][std::max(number_of_panes,3)-2][1].asDouble();
        Real64 eog_c = root["Spacers"][spacer_type]["Coefficients"][std::max(number_of_panes,3)-2][2].asDouble();
        u_eog = eog_a + eog_b*u_cog + eog_c*pow_2(u_cog);
      }

      frame_edge_ratio = u_eog/u_cog;

      // Set frame and divider properties
      if ( has_frame )
      {
        new_frame_divider.Name = construction_name + ":FRAME";
        new_frame_divider.FrameWidth = frame_width;
        new_frame_divider.FrameProjectionOut = 0.0;
        new_frame_divider.FrameProjectionIn = 0.0;
        new_frame_divider.FrameConductance = frame_conductance;
        new_frame_divider.FrEdgeToCenterGlCondRatio = frame_edge_ratio;
        new_frame_divider.FrameSolAbsorp = frame_solar_absorptivity;
        new_frame_divider.FrameVisAbsorp = frame_visible_absorptivity;
        new_frame_divider.FrameEmis = frame_IR_emissivity;
        new_frame_divider.FrameEdgeWidth = edge_width; // 2.5 in
        new_frame_divider.DividerType = DividedLite;
        new_frame_divider.DividerWidth = divider_width;
        new_frame_divider.HorDividers = num_horizontal_dividers;
        new_frame_divider.VertDividers = num_vertical_dividers;
        new_frame_divider.DividerProjectionOut = 0.0;
        new_frame_divider.DividerProjectionIn = 0.0;
        new_frame_divider.DividerConductance = frame_conductance;
        new_frame_divider.DivEdgeToCenterGlCondRatio = frame_edge_ratio;
        new_frame_divider.DividerSolAbsorp = frame_solar_absorptivity;
        new_frame_divider.DividerVisAbsorp = frame_visible_absorptivity;
        new_frame_divider.DividerEmis = frame_IR_emissivity;
        new_frame_divider.DividerEdgeWidth = edge_width; // 2.5 in

        frame_area = fenestration_area - glazing_area;
        SurfaceWindow( 1 ).FrameArea = frame_area;
        SurfaceWindow( 1 ).DividerArea = divider_width*(num_horizontal_dividers*glazing_width + num_vertical_dividers*glazing_height - num_horizontal_dividers*num_vertical_dividers*divider_width);
        Surface( 1 ).Area -= SurfaceWindow( 1 ).DividerArea;
        SurfaceWindow( 1 ).GlazedFrac = Surface( 1 ).Area / ( Surface( 1 ).Area + SurfaceWindow( 1 ).DividerArea );

        FrameDivider( 1 ) = new_frame_divider;

        Surface( 1 ).FrameDivider = 1;
      }

      // Calculate total U-factor
      calc_window_performance(in_air_temp, out_air_temp, wind_speed, solar_incident);

      u_factor = -WinHeatGain(1)/(Surface( 1 ).Area*(in_air_temp - out_air_temp));

      // Set up SHGC conditions
      in_air_temp = s_indoor_temp;
      out_air_temp = s_outdoor_temp;
      wind_speed = s_wind_speed;
      solar_incident = s_solar;

      // Calculate SHGC
      calc_window_performance(in_air_temp, out_air_temp, wind_speed, solar_incident);

      Real64 q_total = WinHeatGain(1);

      // NFRC 201-2014 Equation 8-7
      Real64 q_U = u_factor*Surface( 1 ).Area*(out_air_temp - in_air_temp);

      // NFRC 201-2014 Equation 8-2
      shgc = (q_total - q_U)/(Surface( 1 ).Area*solar_incident);

      Real64 non_opaque_area_fraction = Surface( 1 ).Area/fenestration_area;
      vt = POLYF(1.0,Construct( 1 ).TransVisBeamCoef( 1 ))*non_opaque_area_fraction;

      // if match not obtained adjust inputs

      // Deallocate construction specific arrays
      AWinSurf.deallocate();
      QRadSWwinAbs.deallocate();
      QRadSWwinAbsLayer.deallocate();

      if (!u_factor_set && !shgc_set) target_matched = true;

      u_factor_diff = std::abs(target_u_factor - u_factor);
      shgc_diff = std::abs(target_shgc - shgc);

      if (u_factor_diff < u_factor_match_tolerance && shgc_diff < optical_match_tolerance)
      {
        target_matched = true;
      }

      // adjust properties
      if (! target_matched)
      {
        target_matched = true; // TODO adjust properties
      }

    } // end matching loop

    ASHRAE1588RP_Flag = false;
    KickOffSimulation = true;

    if ( ashrae1588_file_name != "" )
    {
      // Write to file
      Json::Value output_1588;
      Json::StyledStreamWriter writer;
      output_1588["Metadata"]["Name"] = construction_name;
      output_1588["Metadata"]["U-factor"] = u_factor;
      output_1588["Metadata"]["Solar Heat Gain Coefficient"] = shgc;
      output_1588["Metadata"]["Fenestration Type"] = fenestration_type;
      output_1588["Metadata"]["Fenstration Area"] = fenestration_area;
      output_1588["Metadata"]["Fenestration Width"] = fenestration_width;
      output_1588["Metadata"]["Fenestration Height"] = fenestration_height;
      output_1588["Metadata"]["1588-RP Matching"]["U-factor Target"] = target_u_factor;
      output_1588["Metadata"]["1588-RP Matching"]["Solar Heat Gain Coefficient Target"] = target_shgc;
      output_1588["Metadata"]["1588-RP Matching"]["U-factor Difference"] = u_factor_diff;
      output_1588["Metadata"]["1588-RP Matching"]["Solar Heat Gain Coefficient Difference"] = shgc_diff;
      output_1588["Metadata"]["Visible Transmittance"] = vt;
      output_1588["Glazing"]["Number of Panes"] = number_of_panes;
      output_1588["Glazing"]["Area"] = glazing_area;

      for ( int MaterNum = 1; MaterNum <= number_of_new_materials; MaterNum += 2 ) {
        int i = (MaterNum-1)/2;
        output_1588["Glazing"]["Panes"][i]["Thickness"] = Material( MaterNum ).Thickness;
        output_1588["Glazing"]["Panes"][i]["Conductivity"] = Material( MaterNum ).Conductivity;
        if (i == 0) {
          output_1588["Glazing"]["Panes"][i]["Coating"] = glazing_coating;
          output_1588["Glazing"]["Panes"][i]["Tint"] = glazing_tint;
        }
        else {
          output_1588["Glazing"]["Panes"][i]["Coating"] = "NONE";
          output_1588["Glazing"]["Panes"][i]["Tint"] = "CLEAR";
        }
        output_1588["Glazing"]["Panes"][i]["Average Solar Transmittance"] = Material( MaterNum ).Trans;
        output_1588["Glazing"]["Panes"][i]["Average Solar Back Side Reflectance"] = Material( MaterNum ).ReflectSolBeamBack;
        output_1588["Glazing"]["Panes"][i]["Average Solar Front Side Reflectance"] = Material( MaterNum ).ReflectSolBeamFront;
        output_1588["Glazing"]["Panes"][i]["Average Visible Transmittance"] = Material( MaterNum ).TransVis;
        output_1588["Glazing"]["Panes"][i]["Average Visible Back Side Reflectance"] = Material( MaterNum ).ReflectVisBeamBack;
        output_1588["Glazing"]["Panes"][i]["Average Visible Front Side Reflectance"] = Material( MaterNum ).ReflectVisBeamFront;
        output_1588["Glazing"]["Panes"][i]["Average Infrared Transmittance"] = Material( MaterNum ).TransThermal;
        output_1588["Glazing"]["Panes"][i]["Average Infrared Back Side Absorptance"] = Material( MaterNum ).AbsorpThermalBack;
        output_1588["Glazing"]["Panes"][i]["Average Infrared Front Side Absorptance"] = Material( MaterNum ).AbsorpThermalFront;
        int spectral_data_index = Material( MaterNum ).GlassSpectralDataPtr;
        for ( int lam = 1; lam <= SpectralData(spectral_data_index).NumOfWavelengths; lam++) {
          output_1588["Glazing"]["Panes"][i]["Wavelengths"][lam-1] = SpectralData(spectral_data_index).WaveLength( lam );
          output_1588["Glazing"]["Panes"][i]["Transmittance"][lam-1] = SpectralData(spectral_data_index).Trans( lam );
          output_1588["Glazing"]["Panes"][i]["Reflectance (Front)"][lam-1] = SpectralData(spectral_data_index).ReflFront( lam );
          output_1588["Glazing"]["Panes"][i]["Reflectance (Back)"][lam-1] = SpectralData(spectral_data_index).ReflBack( lam );
        }
      }

      int num_gases = root["Gases"][gas_type].size();

      for ( int MaterNum = 2; MaterNum <= number_of_new_materials; MaterNum += 2 ) {
        int i = (MaterNum-2)/2;
        output_1588["Glazing"]["Gaps"][i]["Primary Gas"] = gas_type;
        output_1588["Glazing"]["Gaps"][i]["Thickness"] = gap_thickness;
        for (int gas = 0; gas < num_gases; gas++) {
          output_1588["Glazing"]["Gaps"][i]["Mixture"][gas] = root["Gases"][gas_type][gas];
        }
      }

      output_1588["Glazing"]["Center-of-Glass U-factor"] = u_cog;
      output_1588["Frame and Divider"]["Frame Width"] = frame_width;
      output_1588["Frame and Divider"]["Frame Conductance"] = frame_conductance;
      output_1588["Frame and Divider"]["Frame Material"] = frame_material;
      output_1588["Frame and Divider"]["Spacer Type"] = spacer_type;
      output_1588["Frame and Divider"]["Edge-of-Glass U-factor"] = u_eog;
      output_1588["Frame and Divider"]["Edge-of-Glass Conductance Ratio"] = frame_edge_ratio;
      output_1588["Frame and Divider"]["Frame Area"] = frame_area;
      output_1588["Frame and Divider"]["Divider Width"] = divider_width;
      output_1588["Frame and Divider"]["Number of Vertical Dividers"] = num_vertical_dividers;
      output_1588["Frame and Divider"]["Number of Horizontal Dividers"] = num_horizontal_dividers;

      std::ofstream output_file(ashrae1588_file_name);

      writer.write(output_file, output_1588);
      output_file.close();
    }

    // deallocate temporary arrays
    remove_dummy_variables();

    // Restore Spectral Data list and copy in new spectral data
    {
      new_spectraldata = SpectralData;
      SpectralData.deallocate();
      TotSpectralData = TotSpectralDataSave;
      SpectralData.allocate( TotSpectralData + num_spectral_datasets );
      SpectralData( {1,TotSpectralData} ) = SpectralDataSave;
      SpectralData( {TotSpectralData + 1, TotSpectralData + num_spectral_datasets}) = new_spectraldata;

      SpectralDataSave.deallocate();

      TotSpectralData += num_spectral_datasets;
    }

    // Restore materials list and copy in new materials
    {
      // Apply dirt factor to outermost layer
      if (dirt_factor == 0.0) // Don't know why this is done, but it happens for all window constructions
        Material[0].GlassTransDirtFactor = 1.0;
      else
        Material[0].GlassTransDirtFactor = dirt_factor;

      for (int i = 1; i <= (int)Material.size_; i++) {
        if ( Material( i ).Group == WindowGlass ) {
          Material( i ).GlassSpectralDataPtr += TotSpectralDataSave;
        }
      }


      new_materials = Material;
      new_nominal_R = NominalR;

      Material.deallocate();
      NominalR.deallocate();

      TotMaterials = TotMaterialsSave;

      Material.allocate( TotMaterials + number_of_new_materials);
      NominalR.allocate( TotMaterials + number_of_new_materials);
      Material( {1,TotMaterials} ) = MaterialSave( {1,TotMaterials} );
      NominalR( {1,TotMaterials} ) = NominalRSave( {1,TotMaterials} );
      Material( {TotMaterials + 1, TotMaterials + number_of_new_materials} ) = new_materials;
      NominalR( {TotMaterials + 1, TotMaterials + number_of_new_materials} ) = new_nominal_R;

      MaterialSave.deallocate();
      NominalRSave.deallocate();
    }

    // Restore frame and divider list and copy in new frame and divider
    FrameDivider.deallocate();

    TotFrameDivider = TotFrameDividerSave;

    if ( has_frame )
    {
      FrameDivider.allocate( TotFrameDivider + 1 );
      FrameDivider( {1,TotFrameDivider} ) = FrameDividerSave;
      FrameDivider( TotFrameDivider + 1 ) = new_frame_divider;
      TotFrameDivider += 1;
    }
    else
    {
      FrameDivider.allocate( TotFrameDivider );
      FrameDivider = FrameDividerSave;
    }

    FrameDividerSave.deallocate();


    // Restore construction list and copy in new construction
    {
      Real64 newU = NominalU( 1 );
      Real64 newR = NominalRforNominalUCalculation( 1 );

      Construct.deallocate();
      NominalRforNominalUCalculation.deallocate();
      NominalU.deallocate();

      TotConstructs = TotConstructsSave;

      Construct.allocate( TotConstructs );
      Construct = ConstructSave;
      NominalRforNominalUCalculation.allocate( TotConstructs );
      NominalU.allocate( TotConstructs );
      NominalRforNominalUCalculation = NominalRforNominalUCalculationSave;
      NominalU = NominalUSave;


      Construct( ConstrNum ) = new_construct;
      // Set new layer references corresponding to new material numbers
      for ( int Layer = 1; Layer <= Construct( ConstrNum ).TotLayers; ++Layer ) {
        Construct( ConstrNum ).LayerPoint( Layer ) = TotMaterials + Layer;
      }

      TotMaterials += number_of_new_materials;

      if ( has_frame )
      {
        Construct( ConstrNum ).W5FrameDivider = TotFrameDivider;
      }

      NominalRforNominalUCalculation( ConstrNum ) = newR;
      NominalU( ConstrNum ) = newU;

      ConstructSave.deallocate();
      NominalRforNominalUCalculationSave.deallocate();
      NominalUSave.deallocate();
    }


  } // ...end of WindowASHRAE1588RP Constructions DO loop

  if ( stand_alone_analysis )
  {
    // Write to console
    std::string Elapsed;
    int Hours; // Elapsed Time Hour Reporting
    int Minutes; // Elapsed Time Minute Reporting
    Real64 Seconds; // Elapsed Time Second Reporting
    gio::Fmt ETimeFmt( "(I2.2,'hr ',I2.2,'min ',F5.2,'sec')" );
    std::string NumWarnings = RoundSigDigits( TotalWarningErrors );
    strip( NumWarnings );
    std::string NumSevere = RoundSigDigits( TotalSevereErrors );
    strip( NumSevere );

    Time_Finish = epElapsedTime();
    if ( Time_Finish < Time_Start ) Time_Finish += 24.0 * 3600.0;
    Elapsed_Time = Time_Finish - Time_Start;
    Hours = Elapsed_Time / 3600.0;
    Elapsed_Time -= Hours * 3600.0;
    Minutes = Elapsed_Time / 60.0;
    Elapsed_Time -= Minutes * 60.0;
    Seconds = Elapsed_Time;
    if ( Seconds < 0.0 ) Seconds = 0.0;
    gio::write( Elapsed, ETimeFmt ) << Hours << Minutes << Seconds;
    gio::write( "(A)" ) << ( "EnergyPlus ASHRAE 1588-RP Window Construction(s) Generated Successfully-- Elapsed Time=" + Elapsed );
    ShowMessage( "EnergyPlus ASHRAE 1588-RP Window Construction(s) Generated Successfully-- " + NumWarnings + " Warning; " + NumSevere + " Severe Errors; Elapsed Time=" + Elapsed );

    CloseOutOpenFiles();
    exit (EXIT_SUCCESS);
  }

}

Json::Value read_1588_database(std::string file_path)
{
  Json::Value root;
  Json::Reader json_reader;
  std::ifstream db(file_path, std::ifstream::binary);
  bool read_successful = json_reader.parse(db, root, false);
  if (!read_successful) {
    ShowSevereError( "WindowASHRAE1588RP: Could not open fenestration database file." );
    ShowFatalError( "Program terminates for preceding conditions." );
  }
  db.close();
  return root;
}

void search_database_keys_for_input(const std::string &construction_name, const std::string &field_name, const std::vector< std::string > &keys, const std::string &input, bool &ErrorsFound)
{
  if (!(std::find(keys.begin(), keys.end(), input) != keys.end())) {
    std::string message = "Construction:WindowASHRAE1588RP=" + construction_name + ", " + field_name + "=" + input + " not found in 1588 database.";
    ErrorsFound = true;
    ShowSevereError( message );
  }
}


void calc_window_performance(Real64 T_in, Real64 T_out, Real64 v_ws, Real64 I_s)
{
  InitGlassOpticalCalculations();

  // Calculate window performance
  Surface( 1 ).OutDryBulbTemp = T_out;
  TempEffBulkAir( 1 ) = T_in;

  SurfaceWindow( 1 ).IRfromParentZone = StefanBoltzmann*std::pow(T_in + KelvinConv,4);

  // initial guess temperatures
  int num_temps = 2 + 2*Construct(1).TotGlassLayers;
  Real64 in_surf_temp = T_in - (1.0/(num_temps-1))*(T_in - T_out);
  Real64 out_surf_temp = T_out + (1.0/(num_temps-1))*(T_in - T_out);

  Real64 h_exterior_f = 4 + v_ws*4;
  Real64 h_exterior;

  BeamSolarRad = I_s;

  if ( I_s > 0.0 ) {
    SunIsUp = true;
  }

  InitSolarHeatGains();
  CalcInteriorSolarDistribution();

  // Calculate heat balance (iteratively solve for surface temperatures)
  Real64 out_surf_temp_prev = out_surf_temp;
  Real64 in_surf_temp_prev = in_surf_temp;

  Real64 out_surf_temp_diff;
  Real64 in_surf_temp_diff;

  int max_iterations = 20;
  Real64 tolerance = 0.1; // deg C

  // Save tilt information for natural convection calculations
  Real64 tilt_save = Surface( 1 ).Tilt;

  for (int i = 0; i < max_iterations; i++) {

    // Use complementary angle for exterior natural convection calculations
    Surface( 1 ).Tilt = 180 - tilt_save;
    Surface( 1 ).CosTilt = cos(Surface( 1 ).Tilt*Pi/180);
    Surface( 1 ).SinTilt = sin(Surface( 1 ).Tilt*Pi/180);
    CalcISO15099WindowIntConvCoeff( 1, out_surf_temp, T_out); // This subroutine sets the global HConvIn( 1 ) variable. We will use it to set the exterior natural convection.
    h_exterior = h_exterior_f + HConvIn( 1 ); // add natural convection

    // revert tilt for interior natural convection calculations
    Surface( 1 ).Tilt = tilt_save;
    Surface( 1 ).CosTilt = cos(tilt_save*Pi/180);
    Surface( 1 ).SinTilt = sin(tilt_save*Pi/180);
    CalcISO15099WindowIntConvCoeff( 1, in_surf_temp, T_in); // This time it's actually being used as intended. HConvIn( 1 ) is referenced from the actual heat balance calculation.

    CalcWindowHeatBalance( 1, h_exterior, in_surf_temp, out_surf_temp );

    out_surf_temp_diff = std::fabs(out_surf_temp - out_surf_temp_prev);
    in_surf_temp_diff = std::fabs(in_surf_temp - in_surf_temp_prev);

    if ( (out_surf_temp_diff < tolerance) && (in_surf_temp_diff < tolerance) ) {
      break;
    }

    out_surf_temp_prev = out_surf_temp;
    in_surf_temp_prev = in_surf_temp;

  }

}

void create_dummy_variables()
{

  // Zone
  Zone.allocate(1);
  NumOfZones = 1;
  Zone( 1 ).SurfaceFirst = 1;
  Zone( 1 ).SurfaceLast = 1;

  MAT.allocate(1);
  ZoneAirHumRatAvg.dimension(1, 0.0);
  ZoneAirHumRat.dimension(1, 0.0);
  DSZone.allocate(1);
  DGZone.allocate(1);
  DBZoneSSG.allocate(1);
  DBZone.allocate(1);
  ZoneTransSolar.allocate(1);
  ZoneTransSolarEnergy.allocate(1);
  ZoneBmSolFrExtWinsRep.allocate(1);
  ZoneDifSolFrExtWinsRep.allocate(1);
  ZoneBmSolFrExtWinsRepEnergy.allocate(1);
  ZoneDifSolFrExtWinsRepEnergy.allocate(1);
  ZoneBmSolFrIntWinsRep.allocate(1);
  ZoneBmSolFrIntWinsRepEnergy.allocate(1);

  // Surface
  Surface.allocate(1);
  TotSurfaces = 1;
  SurfaceWindow.allocate(1);
  TotWindows = 1;
  Surface( 1 ).Class = SurfaceClass_Window;
  Surface( 1 ).HeatTransSurf = true;
  // Skip base surface stuff?
  Surface( 1 ).BaseSurf = 1; // It's own base surface?
  Surface( 1 ).ExtBoundCond = 0;
  Surface( 1 ).ExtSolar = true;
  Surface( 1 ).ExtWind = true;
  Surface( 1 ).Zone = 1;
  Surface( 1 ).TAirRef = AdjacentAirTemp;

  SurfaceWindow( 1 ).ShadingFlag = -1;
  SurfaceWindow( 1 ).StormWinFlag = -1;

  HConvIn.allocate(1);
  TempEffBulkAir.allocate(1);
  QHTRadSysSurf.dimension(1, 0.0);
  QHWBaseboardSurf.dimension(1, 0.0);
  QSteamBaseboardSurf.dimension(1, 0.0);
  QElecBaseboardSurf.dimension(1, 0.0);
  CosIncAng.allocate( 1, 1, 1 );
  SunlitFrac.allocate(1, 1, 1);
  AOSurf.allocate(1);
  SunlitFracWithoutReveal.allocate(1,1,1);
  QRadThermInAbs.dimension(1, 0.0);
  AirSkyRadSplit.allocate(1);
  QRadSWOutIncident.allocate(1);
  WinHeatGain.allocate(1);
  WinTransSolar.allocate(1);
  WinGainConvGlazToZoneRep.allocate(1);
  WinGainIRGlazToZoneRep.allocate(1);
  WinGapConvHtFlowRep.allocate(1);
  WinGapConvHtFlowRepEnergy.allocate(1);
  QS.dimension(1, 0.0);
  WinLossSWZoneToOutWinRep.allocate(1);
  WinSysSolTransmittance.allocate(1);
  WinSysSolAbsorptance.allocate(1);
  WinSysSolReflectance.allocate(1);
  InsideGlassCondensationFlag.allocate(1);
  QdotConvOutRep.allocate(1);
  QdotConvOutRepPerArea.allocate(1);
  QConvOutReport.allocate(1);
  QdotRadOutRep.allocate(1);
  QdotRadOutRepPerArea.allocate(1);
  QRadOutReport.allocate(1);
  AISurf.allocate(1);
  AOSurf.allocate(1);
  ISABSF.dimension(1, 0.0);
  BmIncInsSurfIntensRep.allocate(1);
  BmIncInsSurfAmountRep.allocate(1);
  BmIncInsSurfAmountRepEnergy.allocate(1);
  WinBmSolar.allocate(1);
  WinDifSolar.allocate(1);
  WinBmSolarEnergy.allocate(1);
  WinDifSolarEnergy.allocate(1);
  WinTransSolarEnergy.allocate(1);
  WinBmBmSolar.allocate(1);
  WinBmDifSolar.allocate(1);
  WinBmBmSolarEnergy.allocate(1);
  WinBmDifSolarEnergy.allocate(1);
  WinDirSolTransAtIncAngle.allocate(1);
  AnisoSkyMult.dimension(1, 0.0); // This may need to change if NFRC adds a diffuse component for SHGC tests
  CosIncidenceAngle.allocate(1);
  QRadSWOutIncidentBeam.allocate(1);
  QRadSWOutIncidentSkyDiffuse.allocate(1);
  QRadSWOutIncidentGndDiffuse.allocate(1);
  QRadSWOutIncBmToDiffReflGnd.allocate(1);
  QRadSWOutIncSkyDiffReflGnd.allocate(1);
  QRadSWwinAbsTot.allocate(1);
  QRadSWwinAbsTotEnergy.allocate(1);
  SWOutAbsTotalReport.allocate(1);
  SWOutAbsTotalReport.allocate(1);
  WinShadingAbsorbedSolar.allocate(1);
  WinGainFrameDividerToZoneRep.allocate(1);
  InsideFrameCondensationFlag.allocate(1);
  InsideDividerCondensationFlag.allocate(1);

  CosIncAng(1,1,1) = 1.0;
  SunlitFrac(1,1,1) = 1.0;
  SunlitFracWithoutReveal(1,1,1) = 1.0;

}

void remove_dummy_variables()
{
  // Zone
  NumOfZones = 0;
  Zone.deallocate();
  MAT.deallocate();
  ZoneAirHumRatAvg.deallocate();
  ZoneAirHumRat.deallocate();
  DSZone.deallocate();
  DGZone.deallocate();
  DBZoneSSG.deallocate();
  DBZone.deallocate();
  ZoneTransSolar.deallocate();
  ZoneTransSolarEnergy.deallocate();
  ZoneBmSolFrExtWinsRep.deallocate();
  ZoneDifSolFrExtWinsRep.deallocate();
  ZoneBmSolFrExtWinsRepEnergy.deallocate();
  ZoneDifSolFrExtWinsRepEnergy.deallocate();
  ZoneBmSolFrIntWinsRep.deallocate();
  ZoneBmSolFrIntWinsRepEnergy.deallocate();

  // Surface
  Surface.deallocate();
  SurfaceWindow.deallocate();
  TempEffBulkAir.deallocate();
  HConvIn.deallocate();
  QHTRadSysSurf.deallocate();
  QHWBaseboardSurf.deallocate();
  QSteamBaseboardSurf.deallocate();
  QElecBaseboardSurf.deallocate();
  CosIncAng.deallocate();
  SunlitFrac.deallocate();
  AOSurf.deallocate();
  SunlitFracWithoutReveal.deallocate();
  QRadThermInAbs.deallocate();
  AirSkyRadSplit.deallocate();
  QRadSWOutIncident.deallocate();
  WinHeatGain.deallocate();
  WinTransSolar.deallocate();
  WinGainConvGlazToZoneRep.deallocate();
  WinGainIRGlazToZoneRep.deallocate();
  WinGapConvHtFlowRep.deallocate();
  WinGapConvHtFlowRepEnergy.deallocate();
  QS.deallocate();
  WinLossSWZoneToOutWinRep.deallocate();
  WinSysSolTransmittance.deallocate();
  WinSysSolAbsorptance.deallocate();
  WinSysSolReflectance.deallocate();
  InsideGlassCondensationFlag.deallocate();
  QdotConvOutRep.deallocate();
  QdotConvOutRepPerArea.deallocate();
  QConvOutReport.deallocate();
  QdotRadOutRep.deallocate();
  QdotRadOutRepPerArea.deallocate();
  QRadOutReport.deallocate();
  AISurf.deallocate();
  AOSurf.deallocate();
  ISABSF.deallocate();
  BmIncInsSurfIntensRep.deallocate();
  BmIncInsSurfAmountRep.deallocate();
  BmIncInsSurfAmountRepEnergy.deallocate();
  WinBmSolar.deallocate();
  WinDifSolar.deallocate();
  WinBmSolarEnergy.deallocate();
  WinDifSolarEnergy.deallocate();
  WinTransSolarEnergy.deallocate();
  WinBmBmSolar.deallocate();
  WinBmDifSolar.deallocate();
  WinBmBmSolarEnergy.deallocate();
  WinBmDifSolarEnergy.deallocate();
  WinDirSolTransAtIncAngle.deallocate();
  AnisoSkyMult.deallocate();
  CosIncidenceAngle.deallocate();
  QRadSWOutIncidentBeam.deallocate();
  QRadSWOutIncidentSkyDiffuse.deallocate();
  QRadSWOutIncidentGndDiffuse.deallocate();
  QRadSWOutIncBmToDiffReflGnd.deallocate();
  QRadSWOutIncSkyDiffReflGnd.deallocate();
  QRadSWwinAbsTot.deallocate();
  QRadSWwinAbsTotEnergy.deallocate();
  SWOutAbsTotalReport.deallocate();
  SWOutAbsTotalReport.deallocate();
  WinShadingAbsorbedSolar.deallocate();
  WinGainFrameDividerToZoneRep.deallocate();
  InsideFrameCondensationFlag.deallocate();
  InsideDividerCondensationFlag.deallocate();

  // Environment
  BeamSolarRad = 0.0;
  SunIsUp = false;

}

} // WindowASHRAE1588RP

} // EnergyPlus
