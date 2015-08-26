// C++ Headers
#include <string>
#include <iostream>
#include <fstream>

// ObjexxFCL Headers
#include <ObjexxFCL/gio.hh>
#include <ObjexxFCL/Fmath.hh>

// GlobalSearch Headers
#include <genetic_algorithm.hh>
#include <objective_function.hh>

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

	bool standAloneAnalysis = true;	// TODO preprocess this for special executable

	// First read location of the 1588 database file
	if ( GetNumObjectsFound( "DatabaseFile:WindowASHRAE1588RP" ) > 1 ) {
		std::cout << "Error: There should be only one DatabaseFile:WindowASHRAE1588RP object." << std::endl;
		ShowFatalError( "There are multiple DatabaseFile:WindowASHRAE1588RP objects." );
	}

	if ( GetNumObjectsFound( "DatabaseFile:WindowASHRAE1588RP" ) < 1 ) {
		std::cout << "Error: There is no DatabaseFile:WindowASHRAE1588RP object. There must be a DatabaseFile:WindowASHRAE1588RP if your input file contains any Construction:WindowASHRAE1588RP objects." << std::endl;
		ShowFatalError( "There is no DatabaseFile:WindowASHRAE1588RP object. There must be a DatabaseFile:WindowASHRAE1588RP if your input file contains any Construction:WindowASHRAE1588RP objects." );
	}

	int databaseNumAlpha;
	int databaseNumNumeric;
	Array1D_string databaseAlphas( 1 ); // Alpha values array
	Array1D< Real64 > databaseNumerics( 0 ); // Numeric values array
	int IOStat; // IO Status when calling get input subroutine

	GetObjectItem( "DatabaseFile:WindowASHRAE1588RP", 1, databaseAlphas, databaseNumAlpha, databaseNumerics, databaseNumNumeric, IOStat, lNumericFieldBlanks, lAlphaFieldBlanks, cAlphaFieldNames, cNumericFieldNames );

	std::string db1588FilePathInput = databaseAlphas( 1 );
	std::string db1588FilePath;
	bool exists;

	CheckForActualFileName( db1588FilePathInput, exists, db1588FilePath );
	if ( ! exists ) {
		ShowSevereError( "WindowASHRAE1588RP: Could not locate the ASHRAE 1588 window database, expecting it as file name=" + db1588FilePathInput );
		ShowContinueError( "Certain run environments require a full path to be included with the file name in the input field." );
		ShowContinueError( "Try again with putting full path and file name in the field." );
		ShowFatalError( "Program terminates due to these conditions." );
	}

	auto database = ASHRAE1588Database(read1588Database(db1588FilePath));

	// Get list of Coatings from the database
	std::vector < std::string > coatingKeys = database.getTraitOptions("Coatings");

	// Get list of Substrates from the database
	std::vector < std::string > substrateKeys = database.getTraitOptions("Substrates");

	// Get list of Fenestration Types from the database
	std::vector < std::string > typeKeys = database.getTraitOptions("Types");

	// Get list of Frames from the database
	std::vector < std::string > frameKeys = database.getTraitOptions("Frames");

	// Get list of Spacers from the database
	std::vector < std::string > spacerKeys = database.getTraitOptions("Spacers");

	// Get list of Gases from the database
	std::vector < std::string > gasKeys = database.getTraitOptions("Gases");



	int ConstructNumAlpha; // Number of construction alpha names being passed
	int ConstructNumNumeric; // dummy variable for properties being passed
	Array1D_string ConstructAlphas( 8 ); // Construction Alpha names defined
	Array1D< Real64 > ConstructNumerics( 8 ); // Temporary array to transfer construction properties
	bool ErrorInName;
	bool IsBlank;

	int TotWinASHRAE1588Constructs = GetNumObjectsFound( "Construction:WindowASHRAE1588RP" ); // Number of window constructions based on ASHRAE 1588RP

	CurrentModuleObject = "Construction:WindowASHRAE1588RP";

	for ( int Loop = 1; Loop <= TotWinASHRAE1588Constructs; ++Loop ) { // Loop through all WindowASHRAE1588RP constructions.

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

		// Read spectral data from database

		// Save Spectral Data
		int TotSpectralDataSave = TotSpectralData;
		Array1D< SpectralDataProperties > SpectralDataSave;
		SpectralDataSave.allocate( TotSpectralData );
		SpectralDataSave = SpectralData;
		SpectralData.deallocate();


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

		// Name
		std::string constructionName = ConstructAlphas( 1 );

		// U-factor
		Real64 uFactorTarget;
		bool uFactorSet;

		if ( lNumericFieldBlanks( 1 ) )
		{
			uFactorSet = false;
		}
		else
		{
			uFactorSet = true;
			uFactorTarget = ConstructNumerics( 1 );
		}

		// SHGC
		Real64 shgcTarget;
		bool shgcSet;

		if ( lNumericFieldBlanks( 2 ) )
		{
			shgcSet = false;
		}
		else
		{
			shgcSet = true;
			shgcTarget = ConstructNumerics( 2 );
		}

		// set locks as appropriate.
		std::vector<Real64> fenestrationTraits;
		std::vector<bool> variableFlags;
		std::vector<bool> continuousFlags;
		std::vector<double> upperBounds;
		std::vector<double> lowerBounds;

		// Fenestration Type
		if ( lAlphaFieldBlanks( 2 ) )
		{
			fenestrationTraits.push_back(database.getTraitIndexWithMaxUtility("Types"));
		}
		else
		{
			searchDatabaseKeysForInput(constructionName, "Fenestration Type", typeKeys, ConstructAlphas( 2 ), ErrorsFound);
			fenestrationTraits.push_back(database.getTraitIndexByName("Types",ConstructAlphas( 2 )));
		}
		variableFlags.push_back(false); // Fenestration Type is not a variable
		continuousFlags.push_back(database.isContinuous("Types"));
		lowerBounds.push_back(database.getTraitLowerBound("Types"));
		upperBounds.push_back(database.getTraitUpperBound("Types"));

		// Number of Panes
		bool numberOfPanesLock;

		if ( lNumericFieldBlanks( 3 ) )
		{
			numberOfPanesLock = false;
		}
		else
		{
			numberOfPanesLock = true;
			std::string panesName = std::to_string(std::lround(ConstructNumerics( 3 )));
			fenestrationTraits.push_back(database.getTraitIndexByName("Panes",panesName));
		}
		if ( ! numberOfPanesLock )
		{
			fenestrationTraits.push_back(database.getTraitIndexWithMaxUtility("Panes"));
		}
		variableFlags.push_back(!numberOfPanesLock);
		continuousFlags.push_back(database.isContinuous("Panes"));
		lowerBounds.push_back(database.getTraitLowerBound("Panes"));
		upperBounds.push_back(database.getTraitUpperBound("Panes"));

		// Glazing Thickness
		bool glazingThicknessLock;
		if ( lNumericFieldBlanks( 4 ) )
		{
			glazingThicknessLock = false;
		}
		else
		{
			glazingThicknessLock = true;
			fenestrationTraits.push_back(ConstructNumerics( 4 ));
		}

		if ( ! glazingThicknessLock )
		{
			fenestrationTraits.push_back(database.getContinuousTraitValueWithMaxUtility("Glazing Thickness"));
		}
		variableFlags.push_back(!glazingThicknessLock);
		continuousFlags.push_back(database.isContinuous("Glazing Thickness"));
		lowerBounds.push_back(database.getTraitLowerBound("Glazing Thickness"));
		upperBounds.push_back(database.getTraitUpperBound("Glazing Thickness"));

		// Glazing Substrate
		bool glazingSubstrateLock;

		if ( lAlphaFieldBlanks( 3 ) )
		{
			glazingSubstrateLock = false;
		}
		else
		{
			glazingSubstrateLock = true;
			searchDatabaseKeysForInput(constructionName, "Glazing Substrate Type", substrateKeys, ConstructAlphas( 3 ), ErrorsFound);
			fenestrationTraits.push_back(database.getTraitIndexByName("Substrates",ConstructAlphas( 3 )));
		}

		if (! glazingSubstrateLock )
		{
			fenestrationTraits.push_back(database.getTraitIndexWithMaxUtility("Substrates"));
		}
		variableFlags.push_back(!glazingSubstrateLock);
		continuousFlags.push_back(database.isContinuous("Substrates"));
		lowerBounds.push_back(database.getTraitLowerBound("Substrates"));
		upperBounds.push_back(database.getTraitUpperBound("Substrates"));

		// Glazing Coating
		bool glazingCoatingLock;

		if ( lAlphaFieldBlanks( 4 ) )
		{
			glazingCoatingLock = false;
		}
		else
		{
			glazingCoatingLock = true;
			searchDatabaseKeysForInput(constructionName, "Glazing Coating Type", coatingKeys, ConstructAlphas( 4 ), ErrorsFound);
			fenestrationTraits.push_back(database.getTraitIndexByName("Coatings",ConstructAlphas( 4 )));
		}

		if (! glazingCoatingLock )
		{
			fenestrationTraits.push_back(database.getTraitIndexWithMaxUtility("Coatings"));
		}
		variableFlags.push_back(!glazingCoatingLock);
		continuousFlags.push_back(database.isContinuous("Coatings"));
		lowerBounds.push_back(database.getTraitLowerBound("Coatings"));
		upperBounds.push_back(database.getTraitUpperBound("Coatings"));

		// Gas Type
		bool gasTypeLock;

		if ( lAlphaFieldBlanks( 5 ) )
		{
			gasTypeLock = false;
		}
		else
		{
			gasTypeLock = true;
			searchDatabaseKeysForInput(constructionName, "Gas Type", gasKeys, ConstructAlphas( 5 ), ErrorsFound);
			fenestrationTraits.push_back(database.getTraitIndexByName("Gases",ConstructAlphas( 5 )));
		}

		if (! gasTypeLock )
		{
			fenestrationTraits.push_back(database.getTraitIndexWithMaxUtility("Gases"));
		}
		variableFlags.push_back(!gasTypeLock);
		continuousFlags.push_back(database.isContinuous("Gases"));
		lowerBounds.push_back(database.getTraitLowerBound("Gases"));
		upperBounds.push_back(database.getTraitUpperBound("Gases"));

		// Gap Thickness
		bool gapThicknessLock;
		if ( lNumericFieldBlanks( 5 ) )
		{
			gapThicknessLock = false;
		}
		else
		{
			gapThicknessLock = true;
			fenestrationTraits.push_back(ConstructNumerics( 5 ));
		}

		if ( ! gapThicknessLock )
		{
			fenestrationTraits.push_back(database.getContinuousTraitValueWithMaxUtility("Gap Thickness"));
		}
		variableFlags.push_back(!gapThicknessLock);
		continuousFlags.push_back(database.isContinuous("Gap Thickness"));
		lowerBounds.push_back(database.getTraitLowerBound("Gap Thickness"));
		upperBounds.push_back(database.getTraitUpperBound("Gap Thickness"));

		// Spacer Material Type
		bool spacerLock;

		if ( lAlphaFieldBlanks( 6 ) )
		{
			spacerLock = false;
		}
		else
		{
			spacerLock = true;
			searchDatabaseKeysForInput(constructionName, "Spacer Material Type", spacerKeys, ConstructAlphas( 6 ), ErrorsFound);
			fenestrationTraits.push_back(database.getTraitIndexByName("Spacers",ConstructAlphas( 6 )));
		}

		if (! spacerLock )
		{
			fenestrationTraits.push_back(database.getTraitIndexWithMaxUtility("Spacers"));
		}
		variableFlags.push_back(!spacerLock);
		continuousFlags.push_back(database.isContinuous("Spacers"));
		lowerBounds.push_back(database.getTraitLowerBound("Spacers"));
		upperBounds.push_back(database.getTraitUpperBound("Spacers"));

		// Frame Material
		bool frameLock;

		if ( lAlphaFieldBlanks( 7 ) )
		{
			frameLock = false;
		}
		else
		{
			frameLock = true;
			searchDatabaseKeysForInput(constructionName, "Frame Material Type", frameKeys, ConstructAlphas( 7 ), ErrorsFound);
			fenestrationTraits.push_back(database.getTraitIndexByName("Frames",ConstructAlphas( 7 )));
		}

		if (! frameLock )
		{
			fenestrationTraits.push_back(database.getTraitIndexWithMaxUtility("Frames"));
		}
		variableFlags.push_back(!frameLock);
		continuousFlags.push_back(database.isContinuous("Frames"));
		lowerBounds.push_back(database.getTraitLowerBound("Frames"));
		upperBounds.push_back(database.getTraitUpperBound("Frames"));

		// Frame Width
		bool frameWidthLock;

		if ( lNumericFieldBlanks( 6 ) )
		{
			frameWidthLock = false;
		}
		else
		{
			frameWidthLock = true;
			fenestrationTraits.push_back(ConstructNumerics( 6 ));
		}

		if (! frameWidthLock )
		{
			fenestrationTraits.push_back(database.getContinuousTraitValueWithMaxUtility("Frame Width"));
		}
		variableFlags.push_back(!frameWidthLock);
		continuousFlags.push_back(database.isContinuous("Frame Width"));
		lowerBounds.push_back(database.getTraitLowerBound("Frame Width"));
		upperBounds.push_back(database.getTraitUpperBound("Frame Width"));

		// Divider Width
		bool dividerWidthLock;

		if ( lNumericFieldBlanks( 7 ) )
		{
			dividerWidthLock = false;
		}
		else
		{
			dividerWidthLock = true;
			fenestrationTraits.push_back(ConstructNumerics( 7 ));
		}

		if (! dividerWidthLock )
		{
			fenestrationTraits.push_back(0.0);
		}
		variableFlags.push_back(false); // hard code divider width as NOT a variable
		continuousFlags.push_back(true);
		lowerBounds.push_back(0.0);
		upperBounds.push_back(std::min(20e-3, fenestrationTraits.at(7)));

		// Dirt Factor
		Real64 dirtFactor;
		if ( lNumericFieldBlanks( 8 ) )
		{
			dirtFactor = 1.0;
		}
		else
		{
			dirtFactor = ConstructNumerics( 8 );
		}

		std::string ashrae1588FileName;
		if ( lAlphaFieldBlanks( 8 ) )
		{
			ashrae1588FileName = "";
		}
		else
		{
			ashrae1588FileName = ConstructAlphas(8);
		}

	// Check for Errors
	if ( ErrorsFound ) {
		std::cout << "Error found in processing ASHRAE 1588 window construction input." << std::endl;
		ShowFatalError( "Error found in processing ASHRAE 1588 window construction input." );
	}

		// Save Frame and Divider objects
		int TotFrameDividerSave = TotFrameDivider;
		Array1D< FrameDividerProperties > FrameDividerSave;

		FrameDividerSave.allocate( TotFrameDivider );
		FrameDividerSave = FrameDivider;

		FrameDivider.deallocate();

		// Allocate temporary arrays

		ASHRAE1588RP_Flag = true;
		KickOffSimulation = false;

		std::size_t populationSize{50};
		std::size_t maximumGenerations{1000};
		Real64 probabilityOfCrossover{0.8};
		Real64 probabilityOfMutation{0.15};

		auto fsObjFn = FenSysObjFunc(constructionName, database, uFactorTarget, shgcTarget);
		auto ga = globalSearch::GeneticAlgorithm(
				&fsObjFn,
				fenestrationTraits,
				continuousFlags,
				variableFlags,
				upperBounds,
				lowerBounds,
				populationSize,
				maximumGenerations,
				probabilityOfCrossover,
				probabilityOfMutation);
		globalSearch::Individual optimalFenestrationTraits = ga.run();
		auto fs = FenestrationSystem(
				constructionName,
				database,
				uFactorTarget,
				shgcTarget,
				optimalFenestrationTraits.genes);
		fs.calculate();

		ASHRAE1588RP_Flag = false;
		KickOffSimulation = true;

		if ( ashrae1588FileName != "" )
		{
			// Write to file
			Json::StyledStreamWriter writer;

			std::ofstream outputFile(ashrae1588FileName);

			writer.write(outputFile, fs.generateOutput());
			outputFile.close();
		}

		// deallocate temporary arrays

		// Restore Spectral Data list and copy in new spectral data
		{
			Array1D< SpectralDataProperties > newSpectralData = SpectralData;
			int numSpectralDatasets = SpectralData.size_;
			SpectralData.deallocate();
			TotSpectralData = TotSpectralDataSave;
			SpectralData.allocate( TotSpectralData + numSpectralDatasets );
			SpectralData( {1,TotSpectralData} ) = SpectralDataSave;
			SpectralData( {TotSpectralData + 1, TotSpectralData + numSpectralDatasets}) = newSpectralData;

			SpectralDataSave.deallocate();

			TotSpectralData += numSpectralDatasets;
		}

		// Restore materials list and copy in new materials
		// Apply dirt factor to outermost layer
		if (dirtFactor == 0.0) // Don't know why this is done, but it happens for all window constructions
			Material[0].GlassTransDirtFactor = 1.0;
		else
			Material[0].GlassTransDirtFactor = dirtFactor;

		for (int i = 1; i <= (int)Material.size_; i++) {
			if ( Material( i ).Group == WindowGlass ) {
				Material( i ).GlassSpectralDataPtr += TotSpectralDataSave;
			}
		}


		Array1D< MaterialProperties > newMaterials = Material;
		Array1D< Real64 > newNominalR = NominalR;
		int numberOfNewMaterials = Material.size_;

		Material.deallocate();
		NominalR.deallocate();

		TotMaterials = TotMaterialsSave;

		Material.allocate( TotMaterials + numberOfNewMaterials);
		NominalR.allocate( TotMaterials + numberOfNewMaterials);
		Material( {1,TotMaterials} ) = MaterialSave( {1,TotMaterials} );
		NominalR( {1,TotMaterials} ) = NominalRSave( {1,TotMaterials} );
		Material( {TotMaterials + 1, TotMaterials + numberOfNewMaterials} ) = newMaterials;
		NominalR( {TotMaterials + 1, TotMaterials + numberOfNewMaterials} ) = newNominalR;

		MaterialSave.deallocate();
		NominalRSave.deallocate();

		// Restore frame and divider list and copy in new frame and divider
		TotFrameDivider = TotFrameDividerSave;

		bool hasFrame = (FrameDivider.size_ > 0);

		if ( hasFrame )
		{
			FrameDividerProperties newFrameDivider = FrameDivider( 1 );
			FrameDivider.allocate( TotFrameDivider + 1 );
			FrameDivider( {1,TotFrameDivider} ) = FrameDividerSave;
			FrameDivider( TotFrameDivider + 1 ) = newFrameDivider;
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
			ConstructionData newConstruction = Construct( 1 );

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


			Construct( ConstrNum ) = newConstruction;
			// Set new layer references corresponding to new material numbers
			for ( int Layer = 1; Layer <= Construct( ConstrNum ).TotLayers; ++Layer ) {
				Construct( ConstrNum ).LayerPoint( Layer ) = TotMaterials + Layer;
			}

			TotMaterials += numberOfNewMaterials;

			if ( hasFrame )
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

	if ( standAloneAnalysis )
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

Json::Value read1588Database(std::string filePath)
{
	Json::Value root;
	Json::Reader jsonReader;
	std::ifstream db(filePath, std::ifstream::binary);
	bool readSuccessful = jsonReader.parse(db, root, false);
	if (!readSuccessful) {
		ShowSevereError( "WindowASHRAE1588RP: Could not open fenestration database file." );
		ShowFatalError( "Program terminates for preceding conditions." );
	}
	db.close();
	return root;
}

void searchDatabaseKeysForInput(const std::string &constructionName, const std::string &fieldName, const std::vector< std::string > &keys, const std::string &input, bool &ErrorsFound)
{
	if (!(std::find(keys.begin(), keys.end(), input) != keys.end())) {
		std::string message = "Construction:WindowASHRAE1588RP=" + constructionName + ", " + fieldName + "=" + input + " not found in 1588 database.";
		ErrorsFound = true;
		ShowSevereError( message );
	}
}


void calcWindowPerformance(Real64 T_in, Real64 T_out, Real64 v_ws, Real64 I_s)
{
	InitGlassOpticalCalculations();

	// Calculate window performance
	Surface( 1 ).OutDryBulbTemp = T_out;
	TempEffBulkAir( 1 ) = T_in;

	SurfaceWindow( 1 ).IRfromParentZone = StefanBoltzmann*std::pow(T_in + KelvinConv,4);

	// initial guess temperatures
	int numTemps = 2 + 2*Construct(1).TotGlassLayers;
	Real64 inSurfTemp = T_in - (1.0/(numTemps-1))*(T_in - T_out);
	Real64 outSurfTemp = T_out + (1.0/(numTemps-1))*(T_in - T_out);

	Real64 h_exterior_f = 4 + v_ws*4;
	Real64 h_exterior;

	BeamSolarRad = I_s;

	if ( I_s > 0.0 ) {
		SunIsUp = true;
	}

	InitSolarHeatGains();
	CalcInteriorSolarDistribution();

	// Calculate heat balance (iteratively solve for surface temperatures)
	Real64 outSurfTempPrev = outSurfTemp;
	Real64 inSurfTempPrev = inSurfTemp;

	Real64 outSurfTempDiff;
	Real64 inSurfTempDiff;

	int maxIterations = 20;
	Real64 tolerance = 0.1; // deg C

	// Save tilt information for natural convection calculations
	Real64 tiltSave = Surface( 1 ).Tilt;

	for (int i = 0; i < maxIterations; i++) {

		// Use complementary angle for exterior natural convection calculations
		Surface( 1 ).Tilt = 180 - tiltSave;
		Surface( 1 ).CosTilt = cos(Surface( 1 ).Tilt*Pi/180);
		Surface( 1 ).SinTilt = sin(Surface( 1 ).Tilt*Pi/180);
		CalcISO15099WindowIntConvCoeff( 1, outSurfTemp, T_out); // This subroutine sets the global HConvIn( 1 ) variable. We will use it to set the exterior natural convection.
		h_exterior = h_exterior_f + HConvIn( 1 ); // add natural convection

		// revert tilt for interior natural convection calculations
		Surface( 1 ).Tilt = tiltSave;
		Surface( 1 ).CosTilt = cos(tiltSave*Pi/180);
		Surface( 1 ).SinTilt = sin(tiltSave*Pi/180);
		CalcISO15099WindowIntConvCoeff( 1, inSurfTemp, T_in); // This time it's actually being used as intended. HConvIn( 1 ) is referenced from the actual heat balance calculation.

		CalcWindowHeatBalance( 1, h_exterior, inSurfTemp, outSurfTemp );

		outSurfTempDiff = std::fabs(outSurfTemp - outSurfTempPrev);
		inSurfTempDiff = std::fabs(inSurfTemp - inSurfTempPrev);

		if ( (outSurfTempDiff < tolerance) && (inSurfTempDiff < tolerance) ) {
			break;
		}

		outSurfTempPrev = outSurfTemp;
		inSurfTempPrev = inSurfTemp;

	}

}

FenestrationSystem::FenestrationSystem(
	const std::string &constructionName,
	const ASHRAE1588Database &db,
	const Real64 &uFactorTarget,
	const Real64 &shgcTarget,
	std::vector <Real64> fenestrationTraits)
	: database(db),
	constructionName(constructionName),
	uFactorTarget(uFactorTarget),
	shgcTarget(shgcTarget)
{
	fenestrationTypeIndex = int(fenestrationTraits[0]);
	numberOfPanesIndex = int(fenestrationTraits[1]);
	glazingThickness = fenestrationTraits[2];
	glazingSubstrateIndex = int(fenestrationTraits[3]);
	glazingCoatingIndex = int(fenestrationTraits[4]);
	gasTypeIndex = int(fenestrationTraits[5]);
	gapThickness = fenestrationTraits[6];
	spacerTypeIndex = int(fenestrationTraits[7]);
	frameTypeIndex = int(fenestrationTraits[8]);
	frameWidth = fenestrationTraits[9];
	dividerWidth = fenestrationTraits[10];

	fenestrationType = database.getTraitNameByIndex("Types",fenestrationTypeIndex);
	numberOfPanes = database.getTraitValueByIndex("Panes",numberOfPanesIndex).asInt();
	glazingSubstrateType = database.getTraitNameByIndex("Substrates",glazingSubstrateIndex);
	glazingCoatingType = database.getTraitNameByIndex("Coatings",glazingCoatingIndex);
	gasType = database.getTraitNameByIndex("Gases",gasTypeIndex);
	spacerType = database.getTraitNameByIndex("Spacers",spacerTypeIndex);
	frameType = database.getTraitNameByIndex("Frames",frameTypeIndex);


}


void
FenestrationSystem::calculate() {

	createDummyVariables();

	Construct( 1 ).Name = constructionName;
	Construct( 1 ).TypeIsWindow = true;

	Surface( 1 ).Name = constructionName + ":Surface";

	bool ErrorsFound = false;
	int spectralDataSize = database.wavelengths.size();

	// Set testing conditions
	Real64 indoorTempUFactor = database.tests["U-factor"]["Indoor Temperature"].asDouble();
	Real64 outdoorTempUFactor = database.tests["U-factor"]["Outdoor Temperature"].asDouble();
	Real64 windSpeedUFactor = database.tests["U-factor"]["Wind Speed"].asDouble();
	Real64 solarUFactor = database.tests["U-factor"]["Solar Incidence"].asDouble();

	Real64 indoorTempSHGC = database.tests["SHGC"]["Indoor Temperature"].asDouble();
	Real64 outdoorTempSHGC = database.tests["SHGC"]["Outdoor Temperature"].asDouble();
	Real64 windSpeedSHGC = database.tests["SHGC"]["Wind Speed"].asDouble();
	Real64 solarSHGC = database.tests["SHGC"]["Solar Incidence"].asDouble();

	Json::Value fenestrationTypeValue = database.getTraitValueByIndex("Types",fenestrationTypeIndex);
	Json::Value substrateValue = database.getTraitValueByIndex("Substrates",glazingSubstrateIndex);
	Json::Value coatingValue = database.getTraitValueByIndex("Coatings",glazingCoatingIndex);
	Json::Value gasValue = database.getTraitValueByIndex("Gases",gasTypeIndex);
	Json::Value spacerTypeValue = database.getTraitValueByIndex("Spacers",spacerTypeIndex);
	Json::Value frameMaterialValue = database.getTraitValueByIndex("Frames",frameTypeIndex);

	Real64 frameSolarAbsorptivity = fenestrationTypeValue["Frame Absorptivity"].asDouble();
	Real64 frameVisibleAbsorptivity = fenestrationTypeValue["Frame Absorptivity"].asDouble();

	Real64 frameUFactor;
	std::string fenestrationCategory = fenestrationTypeValue["Category"].asString();
	if ( spacerTypeValue["Metal"].asBool() ) {
		frameUFactor = frameMaterialValue["Metal Spacer"][fenestrationCategory][std::min(numberOfPanes,3)-1].asDouble();
	}
	else {
		frameUFactor = frameMaterialValue["Non-Metal Spacer"][fenestrationCategory][std::min(numberOfPanes,3)-1].asDouble();
	}

	Real64 assumed_h_o = 30;	// W/m2-K
	Real64 assumed_h_i = 8;	// W/m2-K

	if ( (1/assumed_h_o + 1/assumed_h_i) >= (1/frameUFactor) ) {
		frameConductance = 9999999;
	}
	else {
		frameConductance = 1/(1/frameUFactor - (1/assumed_h_o + 1/assumed_h_i));
	}

	fenestrationWidth = fenestrationTypeValue["Width"].asDouble();
	fenestrationHeight = fenestrationTypeValue["Height"].asDouble();
	Real64 tilt = fenestrationTypeValue["Tilt"].asDouble()*Pi/180.0;

	bool hasFrame;

	if ( frameWidth > 0.0 )
	{
		hasFrame = true;
	}
	else
	{
		hasFrame = false;
	}

	fenestrationArea = fenestrationWidth*fenestrationHeight;
	Real64 glazingWidth = fenestrationWidth - 2.0*frameWidth;
	Real64 glazingHeight = fenestrationHeight - 2.0*frameWidth;
	glazingArea = glazingWidth*glazingHeight;

	Real64 maxDividerSpacing = 0.3; // NFRC 100-2014 4.2.2 (B)
	Real64 frameIREmissivity = 0.9;
	Real64 edgeWidth = 0.06355;

	if ( hasFrame && dividerWidth > 0.0)
	{
		numHorizontalDividers = ceil(glazingHeight/maxDividerSpacing);
		numVerticalDividers = ceil(glazingWidth/maxDividerSpacing);
	}
	else
	{
		numHorizontalDividers = 0;
		numVerticalDividers = 0;
	}

	Surface( 1 ).Height = glazingHeight;
	Surface( 1 ).Width = glazingWidth;
	Surface( 1 ).Area = glazingArea;
	Surface( 1 ).Tilt = tilt*180/Pi;
	Surface( 1 ).CosTilt = cos(tilt);
	Surface( 1 ).SinTilt = sin(tilt);
	Surface( 1 ).ViewFactorSky = 0.5 * ( 1.0 + Surface( 1 ).CosTilt );
	Surface( 1 ).ViewFactorGround = 0.5 * ( 1.0 - Surface( 1 ).CosTilt );
	Surface( 1 ).ViewFactorSkyIR = Surface( 1 ).ViewFactorSky;
	Surface( 1 ).ViewFactorGroundIR = Surface( 1 ).ViewFactorGround;
	AirSkyRadSplit( 1 ) = std::sqrt( 0.5 * ( 1.0 + Surface( 1 ).CosTilt ) );

	int numberOfGaps = numberOfPanes - 1;
	int numberOfNewMaterials = numberOfPanes + numberOfGaps;

	// Construction specific allocations
	AWinSurf.allocate(numberOfPanes, 1);
	QRadSWwinAbs.allocate(numberOfPanes, 1);
	QRadSWwinAbsLayer.allocate(numberOfPanes, 1);

	// Create New Spectral Data objects
	int numSpectralDatasets = std::min(numberOfPanes,2);
	if ( SpectralData.size_ != (unsigned)numSpectralDatasets ) {
		SpectralData.allocate(numSpectralDatasets);
		TotSpectralData = numSpectralDatasets;
	}

	SpectralData( 1 ).Name = constructionName + ":SPECTRALDATA1";
	SpectralData( 1 ).NumOfWavelengths = spectralDataSize;

	SpectralData( 1 ).WaveLength.deallocate( );
	SpectralData( 1 ).Trans.deallocate( );
	SpectralData( 1 ).ReflFront.deallocate( );
	SpectralData( 1 ).ReflBack.deallocate( );

	SpectralData( 1 ).WaveLength.allocate( spectralDataSize ); // Wavelength (microns)
	SpectralData( 1 ).Trans.allocate( spectralDataSize ); // Transmittance at normal incidence
	SpectralData( 1 ).ReflFront.allocate( spectralDataSize ); // Front reflectance at normal incidence
	SpectralData( 1 ).ReflBack.allocate( spectralDataSize ); // Back reflectance at normal incidence

	for ( int i = 1; i <= spectralDataSize; i++ ) {
		SpectralData( 1 ).WaveLength( i ) = database.wavelengths[i-1]; // Wavelengths

		// Construct spectral data from component properties
		Real64 tau_s = substrateValue["tau_s"][i-1].asDouble();
		tau_s = std::pow(tau_s,glazingThickness/substrateValue["Thickness"].asDouble());
		Real64 r_s = substrateValue["r_s"][i-1].asDouble();
		Real64 t_c, rf_c, rb_c;

		if (database.getTraitNameByIndex("Coatings",glazingCoatingIndex) != "NONE") {
			t_c = coatingValue["t_c"][i-1].asDouble();
			rf_c = coatingValue["rf_c"][i-1].asDouble();
			rb_c = coatingValue["rb_c"][i-1].asDouble();
		}
		else {
			t_c = 1 - r_s;
			rf_c = r_s;
			rb_c = r_s;
		}

		SpectralData( 1 ).Trans( i ) = ((1-r_s)*t_c*tau_s)/(1 - r_s*rf_c*tau_s*tau_s);
		// Following is needed since angular calculation in subr TransAndReflAtPhi
		// fails for Trans = 0.0
		if ( SpectralData( 1 ).Trans( i ) < 0.001 ) {
			SpectralData( 1 ).Trans( i ) = 0.001;
		}
		SpectralData( 1 ).ReflFront( i ) = r_s + (pow_2(1 - r_s)*rf_c*tau_s*tau_s)/(1 - r_s*rf_c*tau_s*tau_s);
		SpectralData( 1 ).ReflBack( i ) = rb_c + (t_c*t_c*r_s*tau_s*tau_s)/(1 - r_s*rf_c*tau_s*tau_s);
	}

	std::string RoutineName = "WindowASHRAE1588RP";

	// Check integrity of the spectral data
	for ( int LamNum = 1; LamNum <= spectralDataSize; ++LamNum ) {
		Real64 Lam = SpectralData( 1 ).WaveLength( LamNum );
		Real64 Tau = SpectralData( 1 ).Trans( LamNum );
		Real64 RhoF = SpectralData( 1 ).ReflFront( LamNum );
		Real64 RhoB = SpectralData( 1 ).ReflBack( LamNum );
		if ( LamNum < spectralDataSize ) {
			if ( SpectralData( 1 ).WaveLength( LamNum + 1 ) <= Lam ) {
				ErrorsFound = true;
				ShowSevereError( RoutineName + CurrentModuleObject + "=\"" + SpectralData( 1 ).Name + "\" invalid set." );
				ShowContinueError( "... Wavelengths not in increasing order. at wavelength#=" + TrimSigDigits( LamNum ) + ", value=[" + TrimSigDigits( Lam, 4 ) + "], next is [" + TrimSigDigits( SpectralData( 1 ).WaveLength( LamNum + 1 ), 4 ) + "]." );
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

	if ( numSpectralDatasets == 2 ) {
		SpectralData( 2 ).Name = constructionName + ":SPECTRALDATA2";
		SpectralData( 2 ).NumOfWavelengths = spectralDataSize;

		SpectralData( 2 ).WaveLength.deallocate( );
		SpectralData( 2 ).Trans.deallocate( );
		SpectralData( 2 ).ReflFront.deallocate( );
		SpectralData( 2 ).ReflBack.deallocate( );

		SpectralData( 2 ).WaveLength.allocate( spectralDataSize ); // Wavelength (microns)
		SpectralData( 2 ).Trans.allocate( spectralDataSize ); // Transmittance at normal incidence
		SpectralData( 2 ).ReflFront.allocate( spectralDataSize ); // Front reflectance at normal incidence
		SpectralData( 2 ).ReflBack.allocate( spectralDataSize ); // Back reflectance at normal incidence

		for ( int i = 1; i <= spectralDataSize; i++ ) {
			SpectralData( 2 ).WaveLength( i ) = SpectralData( 1 ).WaveLength( i ); // Wavelengths

			std::string innerSubstrate = "CLEAR";

			if (database.getTraitNameByIndex("Substrate",glazingSubstrateIndex) == "LOWIRON") {
				innerSubstrate = "LOWIRON";
			}

			Json::Value innerSubstrateValue = database.getTraitValueByName("Substrates",innerSubstrate);

			// Construct spectral data from component properties
			Real64 tau_s = innerSubstrateValue["tau_s"][i-1].asDouble();
			tau_s = std::pow(tau_s,glazingThickness/innerSubstrateValue["Thickness"].asDouble());
			Real64 r_s = innerSubstrateValue["r_s"][i-1].asDouble();
			Real64 t_c = 1 - r_s;
			Real64 rf_c = r_s;
			Real64 rb_c = r_s;

			SpectralData( 2 ).Trans( i ) = ((1-r_s)*t_c*tau_s)/(1 - r_s*rf_c*tau_s*tau_s);
			// Following is needed since angular calculation in subr TransAndReflAtPhi
			// fails for Trans = 0.0
			if ( SpectralData( 2 ).Trans( i ) < 0.001 ) {
				SpectralData( 2 ).Trans( i ) = 0.001;
			}
			SpectralData( 2 ).ReflFront( i ) = r_s + (pow_2(1 - r_s)*rf_c*tau_s*tau_s)/(1 - r_s*rf_c*tau_s*tau_s);
			SpectralData( 2 ).ReflBack( i ) = rb_c + (t_c*t_c*r_s*tau_s*tau_s)/(1 - r_s*rf_c*tau_s*tau_s);
		}

		// Check integrity of the spectral data
		for ( int LamNum = 1; LamNum <= spectralDataSize; ++LamNum ) {
			Real64 Lam = SpectralData( 2 ).WaveLength( LamNum );
			Real64 Tau = SpectralData( 2 ).Trans( LamNum );
			Real64 RhoF = SpectralData( 2 ).ReflFront( LamNum );
			Real64 RhoB = SpectralData( 2 ).ReflBack( LamNum );
			if ( LamNum < spectralDataSize ) {
				if ( SpectralData( 2 ).WaveLength( LamNum + 1 ) <= Lam ) {
					ErrorsFound = true;
					ShowSevereError( RoutineName + CurrentModuleObject + "=\"" + SpectralData( 2 ).Name + "\" invalid set." );
					ShowContinueError( "... Wavelengths not in increasing order. at wavelength#=" + TrimSigDigits( LamNum ) + ", value=[" + TrimSigDigits( Lam, 4 ) + "], next is [" + TrimSigDigits( SpectralData( 2 ).WaveLength( LamNum + 1 ), 4 ) + "]." );
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
	if ( Material.size_ != (unsigned)numberOfNewMaterials ) {
		Material.allocate( numberOfNewMaterials );
		NominalR.allocate( numberOfNewMaterials );
		TotMaterials = numberOfNewMaterials;
	}


	// Define material properties for glazings
	for ( int MaterNum = 1; MaterNum <= numberOfNewMaterials; MaterNum += 2 )
	{
		std::string coating, substrate;
		if (MaterNum == 1) {
			coating = database.getTraitNameByIndex("Coatings",glazingCoatingIndex);
			Material( MaterNum ).GlassSpectralDataPtr = 1;
		}
		else {
			coating = "NONE";
			Material( MaterNum ).GlassSpectralDataPtr = 2;
		}

		Json::Value paneCoatingValue = database.getTraitValueByName("Coatings",coating);

		Material( MaterNum ).Group = WindowGlass;
		Material( MaterNum ).Name = constructionName + ":GLAZING" + std::to_string(MaterNum);
		Material( MaterNum ).Roughness = VerySmooth;
		Material( MaterNum ).ROnly = true;
		Material( MaterNum ).Thickness = glazingThickness;
		Material( MaterNum ).TransThermal = 0.0;
		Material( MaterNum ).AbsorpThermalFront = paneCoatingValue["Emissivity (Front)"].asDouble();
		Material( MaterNum ).AbsorpThermalBack = paneCoatingValue["Emissivity (Back)"].asDouble();
		Material( MaterNum ).Conductivity = 1.0;
		Material( MaterNum ).GlassTransDirtFactor = 1.0;	// Hold at unity to find match and then apply to outside layer
		Material( MaterNum ).YoungModulus = 7.2e10;
		Material( MaterNum ).PoissonsRatio = 0.22;
		Material( MaterNum ).AbsorpThermal = Material( MaterNum ).AbsorpThermalBack;
		Material( MaterNum ).SolarDiffusing = false;

		NominalR( MaterNum ) = Material( MaterNum ).Thickness / Material( MaterNum ).Conductivity;
		Material( MaterNum ).Resistance = NominalR( MaterNum );

	}

	int numGases = gasValue.size();

	// Define material properties for gaps
	for ( int MaterNum = 2; MaterNum <= numberOfNewMaterials; MaterNum += 2 )
	{
		Material( MaterNum ).Group = WindowGasMixture;
		Material( MaterNum ).Name = constructionName + ":GAP" + std::to_string(MaterNum);
		Material( MaterNum ).Roughness = MediumRough;
		Material( MaterNum ).ROnly = true;
		Material( MaterNum ).Thickness = gapThickness;
		Material( MaterNum ).NumberOfGasesInMixture = numGases;

		for ( int gas = 1; gas <= numGases; gas++)
		{
			Material( MaterNum ).GasType( gas ) = 0;
			Material( MaterNum ).GasFract( gas ) = gasValue[gas-1]["Fraction"].asDouble();
			Material( MaterNum ).GasWght( gas ) = gasValue[gas-1]["Molecular Weight"].asDouble();
			Material( MaterNum ).GasSpecHeatRatio( gas ) = gasValue[gas-1]["Specific Heat Ratio"].asDouble();
			for ( int ICoeff = 1; ICoeff <= 3; ++ICoeff ) {
				Material( MaterNum ).GasCon( ICoeff, gas ) = gasValue[gas-1]["Conductivity"][ICoeff-1].asDouble();
				Material( MaterNum ).GasVis( ICoeff, gas ) = gasValue[gas-1]["Viscosity"][ICoeff-1].asDouble();
				Material( MaterNum ).GasCp( ICoeff, gas ) = gasValue[gas-1]["Specific Heat"][ICoeff-1].asDouble();
			}
		}

		Real64 DenomRGas = ( Material( MaterNum ).GasCon( 1, 1 ) + Material( MaterNum ).GasCon( 1, 2 ) * 300.0 + Material( MaterNum ).GasCon( 1, 3 ) * 90000.0 );
		NominalR( MaterNum ) = Material( MaterNum ).Thickness / DenomRGas;

	}

	Construct( 1 ).TotLayers = numberOfNewMaterials;

	for ( int Layer = 1; Layer <= numberOfNewMaterials; ++Layer ) {
		Construct( 1 ).LayerPoint( Layer ) = Layer;
	}

	NominalRforNominalUCalculation( 1 ) = 0.0;
	for ( int Layer = 1; Layer <= Construct( 1 ).TotLayers; ++Layer ) {
		NominalRforNominalUCalculation( 1 ) += NominalR( Construct( 1 ).LayerPoint( Layer ) );
	}

	CheckAndSetConstructionProperties( 1, ErrorsFound );

	Surface( 1 ).Construction = 1; // This is the only construction available to the dummy surface. The actual surface will reference the real construction.
	Surface( 1 ).FrameDivider = 0; // Set temporarily until after Center-of-Glass U-factor is calculated

	// Set up U-factor conditions
	Real64 indoorTemp = indoorTempUFactor;
	Real64 outdoorTemp = outdoorTempUFactor;
	Real64 windSpeed = windSpeedUFactor;
	Real64 solarIncident = solarUFactor;

	// Calculate Center-of-Glass U-factor (without Frame)
	calcWindowPerformance(indoorTempUFactor, outdoorTempUFactor, windSpeedUFactor, solarUFactor);

	uCOG = -WinHeatGain(1)/(Surface( 1 ).Area*(indoorTemp - outdoorTemp));

	if ( numberOfPanes == 1)
		uEOG = uCOG;
	else {
		Real64 eog_a = spacerTypeValue["Coefficients"][std::min(numberOfPanes,3)-2][0].asDouble();
		Real64 eog_b = spacerTypeValue["Coefficients"][std::min(numberOfPanes,3)-2][1].asDouble();
		Real64 eog_c = spacerTypeValue["Coefficients"][std::min(numberOfPanes,3)-2][2].asDouble();
		uEOG = eog_a + eog_b*uCOG + eog_c*pow_2(uCOG);
	}

	frameEdgeRatio = uEOG/uCOG;

	// Set frame and divider properties
	if ( hasFrame )
	{
		FrameDivider.allocate(1);
		TotFrameDivider = 1;

		FrameDivider( 1 ).Name = constructionName + ":FRAME";
		FrameDivider( 1 ).FrameWidth = frameWidth;
		FrameDivider( 1 ).FrameProjectionOut = 0.0;
		FrameDivider( 1 ).FrameProjectionIn = 0.0;
		FrameDivider( 1 ).FrameConductance = frameConductance;
		FrameDivider( 1 ).FrEdgeToCenterGlCondRatio = frameEdgeRatio;
		FrameDivider( 1 ).FrameSolAbsorp = frameSolarAbsorptivity;
		FrameDivider( 1 ).FrameVisAbsorp = frameVisibleAbsorptivity;
		FrameDivider( 1 ).FrameEmis = frameIREmissivity;
		FrameDivider( 1 ).FrameEdgeWidth = edgeWidth; // 2.5 in
		FrameDivider( 1 ).DividerType = DividedLite;
		FrameDivider( 1 ).DividerWidth = dividerWidth;
		FrameDivider( 1 ).HorDividers = numHorizontalDividers;
		FrameDivider( 1 ).VertDividers = numVerticalDividers;
		FrameDivider( 1 ).DividerProjectionOut = 0.0;
		FrameDivider( 1 ).DividerProjectionIn = 0.0;
		FrameDivider( 1 ).DividerConductance = frameConductance;
		FrameDivider( 1 ).DivEdgeToCenterGlCondRatio = frameEdgeRatio;
		FrameDivider( 1 ).DividerSolAbsorp = frameSolarAbsorptivity;
		FrameDivider( 1 ).DividerVisAbsorp = frameVisibleAbsorptivity;
		FrameDivider( 1 ).DividerEmis = frameIREmissivity;
		FrameDivider( 1 ).DividerEdgeWidth = edgeWidth; // 2.5 in

		frameArea = fenestrationArea - glazingArea;
		SurfaceWindow( 1 ).FrameArea = frameArea;
		SurfaceWindow( 1 ).DividerArea = dividerWidth*(numHorizontalDividers*glazingWidth + numVerticalDividers*glazingHeight - numHorizontalDividers*numVerticalDividers*dividerWidth);
		Surface( 1 ).Area -= SurfaceWindow( 1 ).DividerArea;
		SurfaceWindow( 1 ).GlazedFrac = Surface( 1 ).Area / ( Surface( 1 ).Area + SurfaceWindow( 1 ).DividerArea );

		Surface( 1 ).FrameDivider = 1;
	}

	// Calculate total U-factor
	calcWindowPerformance(indoorTemp, outdoorTemp, windSpeed, solarIncident);

	uFactor = -WinHeatGain(1)/(fenestrationArea*(indoorTemp - outdoorTemp));

	// Set up SHGC conditions
	indoorTemp = indoorTempSHGC;
	outdoorTemp = outdoorTempSHGC;
	windSpeed = windSpeedSHGC;
	solarIncident = solarSHGC;

	// Calculate SHGC
	calcWindowPerformance(indoorTemp, outdoorTemp, windSpeed, solarIncident);

	Real64 qTotal = WinHeatGain(1);

	// NFRC 201-2014 Equation 8-7
	Real64 qU = uFactor*fenestrationArea*(outdoorTemp - indoorTemp);

	// NFRC 201-2014 Equation 8-2
	shgc = (qTotal - qU)/(fenestrationArea*solarIncident);

	Real64 nonOpaqueAreaFraction = Surface( 1 ).Area/fenestrationArea;
	visibleTransmittance = POLYF(1.0,Construct( 1 ).TransVisBeamCoef( 1 ))*nonOpaqueAreaFraction;

	// if match not obtained adjust inputs

	// Deallocate construction specific arrays
	AWinSurf.deallocate();
	QRadSWwinAbs.deallocate();
	QRadSWwinAbsLayer.deallocate();

	uFactorDiff = uFactorTarget - uFactor;
	shgcDiff = shgcTarget - shgc;

	error = std::sqrt(pow_2(uFactorDiff) + pow_2(shgcDiff));

	Real64 uFactorMatchTolerance = 0.05; // Precision of NFRC reporting TODO expose?
	Real64 opticalMatchTolerance = 0.01; // Precision of NFRC reporting TODO expose?

	matched = (std::abs(uFactorDiff) < uFactorMatchTolerance) && (std::abs(shgcDiff) < opticalMatchTolerance);
	removeDummyVariables();

} // calculate

Json::Value
FenestrationSystem::generateOutput() {
	Json::Value output1588;
	output1588["Metadata"]["Name"] = constructionName;
	output1588["Metadata"]["U-factor"] = uFactor;
	output1588["Metadata"]["Solar Heat Gain Coefficient"] = shgc;
	output1588["Metadata"]["Fenestration Type"] = fenestrationType;
	output1588["Metadata"]["Fenstration Area"] = fenestrationArea;
	output1588["Metadata"]["Fenestration Width"] = fenestrationWidth;
	output1588["Metadata"]["Fenestration Height"] = fenestrationHeight;
	output1588["Metadata"]["1588-RP Matching"]["U-factor Target"] = uFactorTarget;
	output1588["Metadata"]["1588-RP Matching"]["Solar Heat Gain Coefficient Target"] = shgcTarget;
	output1588["Metadata"]["1588-RP Matching"]["U-factor Difference"] = uFactorDiff;
	output1588["Metadata"]["1588-RP Matching"]["Solar Heat Gain Coefficient Difference"] = shgcDiff;
	output1588["Metadata"]["1588-RP Matching"]["Matched?"] = matched;
	output1588["Metadata"]["Visible Transmittance"] = visibleTransmittance;
	output1588["Glazing"]["Number of Panes"] = numberOfPanes;
	output1588["Glazing"]["Area"] = glazingArea;

	// Since wavelengths are the same for each spectral dataset, they only need
	// to appear once in the output file.
	for ( int lam = 1; lam <= int(database.wavelengths.size()); lam++) {
		output1588["Glazing"]["Wavelengths"][lam-1] = SpectralData( 1 ).WaveLength( lam );
	}

	for ( int MaterNum = 1; MaterNum <= numberOfPanes*2 - 1; MaterNum += 2 ) {
		int i = (MaterNum-1)/2;
		output1588["Glazing"]["Panes"][i]["Thickness"] = Material( MaterNum ).Thickness;
		output1588["Glazing"]["Panes"][i]["Conductivity"] = Material( MaterNum ).Conductivity;
		if (i == 0) {
			output1588["Glazing"]["Panes"][i]["Coating"] = glazingCoatingType;
			output1588["Glazing"]["Panes"][i]["Substrate"] = glazingSubstrateType;
		}
		else {
			output1588["Glazing"]["Panes"][i]["Coating"] = "NONE";
			output1588["Glazing"]["Panes"][i]["Substrate"] = "CLEAR";
		}
		output1588["Glazing"]["Panes"][i]["Average Infrared Transmittance"] = Material( MaterNum ).TransThermal;
		output1588["Glazing"]["Panes"][i]["Average Infrared Back Side Absorptance"] = Material( MaterNum ).AbsorpThermalBack;
		output1588["Glazing"]["Panes"][i]["Average Infrared Front Side Absorptance"] = Material( MaterNum ).AbsorpThermalFront;

		int spectralDataIndex = Material( MaterNum ).GlassSpectralDataPtr;

		Real64 TransSolUp = 0.0,
					 TransVisUp = 0.0,
					 RefFrontSolUp = 0.0,
					 RefFrontVisUp = 0.0,
					 RefBackSolUp = 0.0,
					 RefBackVisUp = 0.0,
					 SolDown = 0.0,
					 VisDown = 0.0;

		for ( int lam = 1; lam <= int(database.wavelengths.size()); lam++) {
			Real64 Trans = SpectralData(spectralDataIndex).Trans( lam );
			Real64 RefFront = SpectralData(spectralDataIndex).ReflFront( lam );
			Real64 RefBack = SpectralData(spectralDataIndex).ReflBack( lam );
			output1588["Glazing"]["Panes"][i]["Transmittance"][lam-1] = Trans;
			output1588["Glazing"]["Panes"][i]["Reflectance (Front)"][lam-1] = RefFront;
			output1588["Glazing"]["Panes"][i]["Reflectance (Back)"][lam-1] = RefBack;

			// Numeric integration of average properties. Follows same method as WindowManager, but uses consistent wavelengths.
			if (lam != int(database.wavelengths.size())) {
				// Spectral Properties
				Real64 TransNext = SpectralData(spectralDataIndex).Trans( lam + 1 );
				Real64 RefFrontNext = SpectralData(spectralDataIndex).ReflFront( lam + 1 );
				Real64 RefBackNext = SpectralData(spectralDataIndex).ReflBack( lam + 1);

				// Wavelengths
				Real64 Wl = SpectralData(spectralDataIndex).WaveLength( lam );
				Real64 WlNext = SpectralData(spectralDataIndex).WaveLength( lam + 1 );

				// Solar Spectrum and Photopic Response
				Real64 SS = database.solarSpectrum[lam - 1];
				Real64 SSNext = database.solarSpectrum[lam];

				Real64 PR = database.photopicResponse[lam - 1];
				Real64 PRNext = database.photopicResponse[lam];

				Real64 eSol = (WlNext - Wl)*0.5*(SS + SSNext);
				Real64 eVis = (WlNext - Wl)*0.5*(SS + SSNext)*0.5*(PR + PRNext);

				TransSolUp += eSol*0.5*(Trans + TransNext);
				RefFrontSolUp += eSol*0.5*(RefFront + RefFrontNext);
				RefBackSolUp += eSol*0.5*(RefBack + RefBackNext);

				TransVisUp += eVis*0.5*(Trans + TransNext);
				RefFrontVisUp += eVis*0.5*(RefFront + RefFrontNext);
				RefBackVisUp += eVis*0.5*(RefBack + RefBackNext);

				SolDown += eSol;
				VisDown += eVis;

			}

		}

		output1588["Glazing"]["Panes"][i]["Average Solar Transmittance"] = TransSolUp/SolDown;
		output1588["Glazing"]["Panes"][i]["Average Solar Front Side Reflectance"] = RefFrontSolUp/SolDown;
		output1588["Glazing"]["Panes"][i]["Average Solar Back Side Reflectance"] = RefBackSolUp/SolDown;
		output1588["Glazing"]["Panes"][i]["Average Visible Transmittance"] = TransVisUp/VisDown;
		output1588["Glazing"]["Panes"][i]["Average Visible Front Side Reflectance"] = RefFrontVisUp/VisDown;
		output1588["Glazing"]["Panes"][i]["Average Visible Back Side Reflectance"] = RefBackVisUp/VisDown;
	}

	for ( int MaterNum = 2; MaterNum <= numberOfPanes*2 - 1; MaterNum += 2 ) {
		int i = (MaterNum-2)/2;
		output1588["Glazing"]["Gaps"][i]["Primary Gas"] = gasType;
		output1588["Glazing"]["Gaps"][i]["Thickness"] = Material( MaterNum ).Thickness;
		for (int gas = 1; gas <= Material( MaterNum ).NumberOfGasesInMixture; gas++) {
			output1588["Glazing"]["Gaps"][i]["Mixture"][gas-1]["Gas"] = database.getTraitValueByIndex("Gases",gasTypeIndex)[gas-1]["Gas"];
			output1588["Glazing"]["Gaps"][i]["Mixture"][gas-1]["Fraction"] = Material( MaterNum ).GasFract( gas );
			output1588["Glazing"]["Gaps"][i]["Mixture"][gas-1]["Molecular Weight"] = Material( MaterNum ).GasWght( gas );
			output1588["Glazing"]["Gaps"][i]["Mixture"][gas-1]["Specific Heat Ratio"] = Material( MaterNum ).GasSpecHeatRatio( gas );
			for ( int ICoeff = 1; ICoeff <= 3; ++ICoeff ) {
				output1588["Glazing"]["Gaps"][i]["Mixture"][gas-1]["Conductivity"][ICoeff-1] = Material( MaterNum ).GasCon( gas, ICoeff );
				output1588["Glazing"]["Gaps"][i]["Mixture"][gas-1]["Viscosity"][ICoeff-1] = Material( MaterNum ).GasVis( gas, ICoeff );
				output1588["Glazing"]["Gaps"][i]["Mixture"][gas-1]["Specific Heat"][ICoeff-1] = Material( MaterNum ).GasCp( gas, ICoeff );
			}
		}
	}

	output1588["Glazing"]["Center-of-Glass U-factor"] = uCOG;
	output1588["Frame and Divider"]["Frame Width"] = frameWidth;
	output1588["Frame and Divider"]["Frame Conductance"] = frameConductance;
	output1588["Frame and Divider"]["Frame Material"] = frameType;
	output1588["Frame and Divider"]["Spacer Type"] = spacerType;
	output1588["Frame and Divider"]["Edge-of-Glass U-factor"] = uEOG;
	output1588["Frame and Divider"]["Edge-of-Glass Conductance Ratio"] = frameEdgeRatio;
	output1588["Frame and Divider"]["Frame Area"] = frameArea;
	output1588["Frame and Divider"]["Divider Width"] = dividerWidth;
	output1588["Frame and Divider"]["Number of Vertical Dividers"] = numVerticalDividers;
	output1588["Frame and Divider"]["Number of Horizontal Dividers"] = numHorizontalDividers;

	return output1588;
}

void createDummyVariables()
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

void removeDummyVariables()
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

ASHRAE1588Database::ASHRAE1588Database(Json::Value db) : db(db) {
	tests = db["Tests"];

	Json::Value &sd = db["Spectral Functions"];
	assert((sd["Wavelengths"].size() == sd["Solar Spectrum"].size()) && (sd["Wavelengths"].size() == sd["Photopic Response"].size()));
	for ( int i = 0; i < int(sd["Wavelengths"].size()); ++i ) {
		wavelengths.push_back(sd["Wavelengths"][i].asDouble());
		solarSpectrum.push_back(sd["Solar Spectrum"][i].asDouble());
		photopicResponse.push_back(sd["Photopic Response"][i].asDouble());
	}
}

std::vector <std::string> ASHRAE1588Database::getTraitOptions(std::string trait) {
	// Only works for discrete traits should add error message if it's not
	std::vector <std::string> traitOptions;
	Json::Value &traitNode = db["Traits"][trait];
	for (Json::ValueIterator itr = traitNode.begin(); itr != traitNode.end(); ++itr) {
		traitOptions.push_back((*itr)["Name"].asString());
	}
	return traitOptions;
}

int ASHRAE1588Database::getTraitIndexByName(std::string trait, std::string name) {
	Json::Value &traitNode = db["Traits"][trait];
	for (struct {Json::ValueIterator itr; int i ;} d = {traitNode.begin(), 0}; d.itr != traitNode.end(); ++d.itr, ++d.i) {
		if ((*d.itr)["Name"].asString() == name) {
			return d.i;
		}
	}
	return -1;
}

std::string ASHRAE1588Database::getTraitNameByIndex(std::string trait, int index) {
	Json::Value &traitNode = db["Traits"][trait];
	return traitNode[index]["Name"].asString();
}

int ASHRAE1588Database::getTraitIndexWithMaxUtility(std::string trait) {
	Json::Value &traitNode = db["Traits"][trait];
	Json::ValueIterator maxItr = std::max_element(traitNode.begin(),traitNode.end(),[](Json::Value a, Json::Value b){return (a["Utility"] < b["Utility"]);});
	return std::distance(traitNode.begin(), maxItr);
}

std::string ASHRAE1588Database::getTraitNameWithMaxUtility(std::string trait) {
	return getTraitNameByIndex(trait, getTraitIndexWithMaxUtility(trait));
}

Json::Value ASHRAE1588Database::getDiscreteTraitValueWithMaxUtility(std::string trait) {
	return getTraitValueByIndex(trait, getTraitIndexWithMaxUtility(trait));
}

Real64 ASHRAE1588Database::getContinuousTraitValueWithMaxUtility(std::string trait) {
	Json::Value &traitNode = db["Traits"][trait];
	if (traitNode["Distribution"].asString() == "Normal") {
		return traitNode["Shift"].asDouble();
	}
	else if (traitNode["Distribution"].asString() == "Log-Normal") {
		Real64 shift = traitNode["Shift"].asDouble();
		Real64 shape = traitNode["Shape"].asDouble();
		return exp(shift-pow_2(shape));
	}
	else {
		return -1;
	}
}

Json::Value ASHRAE1588Database::getTraitValueByName(std::string trait, std::string name) {
	Json::Value &traitNode = db["Traits"][trait];
	int index = getTraitIndexByName(trait, name);
	if (index < 0) {
		std::cerr << "ERROR: " + name + " not found in list of " + trait + "." << std::endl;
	}
	assert(index >= 0);
	return traitNode[index]["Value"];
}

Json::Value ASHRAE1588Database::getTraitValueByIndex(std::string trait, int index) {
	return db["Traits"][trait][index]["Value"];
}

bool
ASHRAE1588Database::isContinuous(std::string trait)
{
	Json::Value &traitNode = db["Traits"][trait];
	return traitNode.isMember("Distribution");
}

Real64
ASHRAE1588Database::getTraitLowerBound(std::string trait)
{
	if (isContinuous(trait)) {
		Json::Value &traitNode = db["Traits"][trait];
		return traitNode["Min"].asDouble();
	}
	return 0.0;
}

Real64
ASHRAE1588Database::getTraitUpperBound(std::string trait)
{
	Json::Value &traitNode = db["Traits"][trait];
	if (isContinuous(trait)) {
		return traitNode["Max"].asDouble();
	}
	return Real64(traitNode.size() - 1);
}

FenSysObjFunc::FenSysObjFunc(
		const std::string& constructionName,
		const ASHRAE1588Database& database,
		const Real64& uFactorTarget,
		const Real64& shgcTarget) :
	constructionName_(constructionName),
	database_(database),
	uFactorTarget_(uFactorTarget),
	shgcTarget_(shgcTarget)
{}

Real64
FenSysObjFunc::call(const std::vector<Real64>& fenestrationTraits)
{
	auto fs = FenestrationSystem(constructionName_, database_, uFactorTarget_, shgcTarget_, fenestrationTraits);
	fs.calculate();
	return -1.0 * fs.error;
}

// void
// setupField(
// 		const int fieldIdx,
// 	 	const bool isAlpha,
// 		const std::string& databaseKey,
// 		const std::string& traitName,
// 		const std::string& constructionName,
// 		const ASHRAE1588Database& database,
// 	 	std::vector<Real64>& fenestrationTraits,
// 		std::vector<bool> variableFlags,
// 		std::vector<bool> continuousFlags
// )
// {
// 	bool lock;
// 	bool continuous = database.isContinuous(traitName);
// 	if (isAlpha) {
// 		if ( lAlphaFieldBlanks( fieldIdx ) ) {
// 			lock = false;
// 			if ( continuous ) {
// 				fenestrationTraits.push_back(database.getContinuousTraitValueWithMaxUtility(traitName));
// 			} else {
// 				fenestrationTraits.push_back(database.getTraitIndexWithMaxUtility(traitName));
// 			}
// 		} else {
// 			lock = true;
// 			searchDatabaseKeysForInput(constructionName, databaseKey, typeKeys, ConstructAlphas( fieldIdx ), ErrorsFound);
// 			fenestrationTraits.push_back(database.getTraitIndexByName(traitName, ConstructAlphas( fieldIdx )));
// 		}
// 	} else {
// 		if ( lNumericFieldBlanks( fieldIdx ) ) {
// 			lock = false;
// 			if ( isContinuous ) {
// 				fenestrationTraits.push_back(database.getContinuousTraitValueWithMaxUtility(traitName));
// 			} else {
// 				fenestrationTraits.push_back(database.getTraitIndexWithMaxUtility(traitName));
// 			}
// 		} else {
// 			lock = true;
// 			if ( continuous ) {
// 				std::string name = std::to_string(std::lround(ConstructNumerics( fieldIdx )));
// 				fenestrationTraits.push_back(database.getTraitIndexByName(traitName, name));
// 			} else {
// 				fenestrationTraits.push_back(ConstructNumerics(fieldIdx));
// 			}
// 		}
// 	}
// 	variableFlags.push_back(!lock);
//	continuousFlags.push_back(database.isContinuous(traitName));
//	lowerBounds.push_back(database.getTraitLowerBound(traitName));
//	upperBounds.push_back(database.getTraitUpperBound(traitName));
// }

} // WindowASHRAE1588RP

} // EnergyPlus
