/**********************************************************************************
	PanTompkins.c
	Author: Hooman Sedghamiz, Jan 2020
	Last Updated: Feb 2021

	-------------------------------------------
	-------------------------------------------
	Description:

	This is an efficient and fixed-point implementation of the well-known Pan and 
	Tompkins beat detector [1]. This code implements an updated version of the 
	published article [2] where the authors published corrections to a few filters
	employed in the paper. Below is a list of stages that the algorithm goes 
	through:

	1. Preprocessing :

		- In preprocessing stage, the signal is passed through 5 filters namely
		low pass (~ <15 Hz), high pass (~ > 5Hz), 5 point derivative, squaring
		and moving average window respectively.

		- The implementation of these filters is very cruicial both in terms
		of end results and computational load. Therefore, the following optimizations
		are performed:

			- No use of convolution: Filters are implemented based on their 
			difference equation making them recursive and minimizing the number 
			of operations. 

			- Additional care has been given to overflow control, while there is
			no garantee that an oveflow would never happen but several gain control
			checkpoints are implemented. Note that still the range and resolution 
			of ADC plays a crucial point in filter outputs.

			- Almost no division nor multiplication is employed and all the operations
			are performed with bit shifting, making the filters extremely efficient.

			- Filters are implemented in both Direct I and Direct II form. 

	2. Signal analysis and peak detection:

		- After signal enhancement, both integrated signal (output of moving average filter)
		and band-passed signal (output of HP filter) are simultanously analyzed for
		beat detection. One of the advantages of this code is the implementation of
		search-back and T-wave discrimination without using any buffers as described 
		below:
			- Search-back: 
				- Once no beat is detected for 166% of the RR mean estimate, a search 
				back needs to be initiated to find the potentially missed beat. In order
				to avoid storing a buffer, the code only stores the tallest peak which
				was previously classified as noise, once there is a need for search-back
				the stored peak is compared against thresholds and therefore no buffer
				is requried.

			- T-wave discrimination:
				- If a beat is close than 360msec to the previous detected beat, a slope
				test to the previous beat is employed to check whether it is a T wave.
				In order to avoid storing a buffer, the code saves the largest peak in 
				the differentiated signal every time a beat is detected. Once a new beat
				is suspicious of being T wave, that stored slope is employed.

	3. Decision making:

		- The algorithms stores 2 seperate thresholds namely one for bandpassed signal and
		one for the integrated. Each of these thresholds is adaptively updated based on
		detecting noise and signal peaks. Once a peak is above threshold in integrated signal
		and also the corresponding peak in BP signal is also above threshold that peak is 
		classified as a beat.
		
	4. Usage:
		- For example usage please see PanTompkinsCMD.c. For use simply include the header 
		and initialize the code (PT_init()). Then simply pass the aquired samples iteratively
		to the main module PT_StateMachine(sample). See below for details.


	References:

	[1] Pan, Jiapu; Tompkins, Willis J. (March 1985). "A Real-Time QRS Detection
	Algorithm". IEEE Transactions on Biomedical Engineering. BME-32 (3): 230�236.

	[2] "A real time QRS detection algorithm (errata corrige), Link:
	https://courses.cs.washington.edu/courses/cse474/18wi/labs/l8/QRSdetection.pdf

	[3] Sedghamiz, Hooman. "Complete Pan Tompkins Implementation ECG QRS detector
	- File Exchange - MATLAB Central" Link:
	https://ww2.mathworks.cn/matlabcentral/fileexchange/45840-complete-pan-tompkins-implementation-ecg-qrs-detector

 **********************************************************************************/


/********************************************************************************
    Headers
 ********************************************************************************/

#include "PanTompkins.h"


/********************************************************************************
    PT Algorithm data and buffers, see PT_init for detailed use of each parameter.
 ********************************************************************************/

static struct PT_struct PT_data;
static struct PT_struct *const PT_dptr = &PT_data;


static int16_t Prev_valBP, 
Prev_Prev_valBP, Best_PeakBP, Prev_valDR, Prev_Prev_valDR, 
Best_PeakDR, Old_PeakDR, Count_SinceRR, RR1_p, RR2_p, 
RR1_sum, RR2_sum, BlankTimeCnt, SBcntI, SB_peakBP, SB_peakDR, 
y_h, st_mean_pkBP;

static uint16_t MV_sum, PEAKI_temp, st_mx_pk, st_mean_pk,
Prev_val, Prev_Prev_val, SB_peakI;


#if (FILTER_FORM == 2)
static int16_t LP_y_new, LP_y_old;
#endif


/**********************************************************************************

    Fuction Name: PT_init


    Parameter:
     Input:   none

     Returns: none

    Description: initializes the PanTompkins (PT) data structure and RR interval,
	and filter Buffers.

 *******************************************************************************/

void PT_init( void )
{
	/**************************************************
	Initialize Pan_Tompkins structure.
	**************************************************/

	memset(&PT_data, 0, sizeof(PT_data));

	PT_dptr->PT_state		= START_UP;

	PT_dptr->Recent_RR_M = PT_dptr->RR_M =  PT1000MS;

	PT_dptr->RR_Low_L		= RR92PERCENT;
	PT_dptr->RR_High_L		= RR116PERCENT;
	PT_dptr->RR_Missed_L	= RR166PERCENT;

	PT_dptr->LP_pointer		= 0;
	PT_dptr->HP_pointer		= 0;
	PT_dptr->MVA_pointer	= 0;

	PT_dptr->HR_State = REGULAR_HR;


	/**************************************************
	Initialize filter buffers
	**************************************************/
	int8_t idex;

	for (idex = 0; idex < LP_BUFFER_SIZE; idex++)
		PT_dptr->LP_buf[idex]		= 0;							//  LP filter buffer
	for (idex = 0; idex < HP_BUFFER_SIZE; idex++)
		PT_dptr->HP_buf[idex]		= 0;							//  HP filter buffer
	for (idex = 0; idex < DR_BUFFER_SIZE; idex++)
		PT_dptr->DR_buf[idex]		= 0;							//  DR filter buffer
	for (idex = 0; idex < MVA_BUFFER_SIZE; idex++)
		PT_dptr->MVA_buf[idex]		= 0;							//  MVA filter buffer
	for (idex = 0; idex < RR_BUFFER_SIZE; idex++) {
		PT_dptr->RR_AVRG1_buf[idex] = 
			PT_dptr->RR_AVRG2_buf[idex] = PT1000MS;					//  Normal	and extreme RR buffers
	}

	/**************************************************
	Initialize all static variables 
	**************************************************/
	Prev_val = Prev_Prev_val = 0;									// Place holders for peak detector in Integrated Sig
	Prev_valBP = Prev_Prev_valBP = Best_PeakBP = 0;					// Place holders for peak detector in BP signal
	Prev_valDR = Prev_Prev_valDR = Best_PeakDR = Old_PeakDR = 0;	// Place holders for peak detector in Derivative signal (Used for T-wave discrimination)
	Count_SinceRR = 0;												// Nr of samples since last qrs peak
	RR1_p = RR2_p = 0;												// Pointers to RR average 1 and 2 resepectively
	MV_sum = 0;														// sum for moving average filter
	RR1_sum = RR2_sum = PT1000MS << 3;								// Sum of RR1 and RR2 buffers
	BlankTimeCnt = 0;												// Counter for blank-time.
	SBcntI = 0;														// For searchback index in Integ Signal
	SB_peakI = 0;													// For searchback in Integ sig
	SB_peakBP = SB_peakDR = 0;										// For searchback peak holders in BP and slope signal
	st_mx_pk = 0;													// Used in learning phase 1 to estimate thresholds
	y_h = 0;														// recusrively used in HP filter

#if (FILTER_FORM == 2)
	LP_y_new = LP_y_old = 0;										// Parameters for DirectForm || LP filter
#endif
}

/**********************************************************************************

	Fuction Name: PT_StateMachine

	Parameter:
	 Input	:	datum		- Most recent sample of ECG from ADC.

	 Returns:	BeatDelay	- If non-zero a qrs has been detected with BeatDelay samples.

	Description: Pan-Tompkins State Machine implementation. This state-machine goes through
	each step of the Pan-Tompkins algorithm, namely, BP filtering, derivative filtering, 
	squaring and integration. After signal enhancement, BP and integrated signal's peaks 
	are analyzed and compared to two adaptive thresholds to determine whether a beat has occured.
	The thresholds themselves is computed adaptively from noise and signal estimation. An efficient
	search-back strategy and T-wave discrimination is also implemented that do not require any 
	buffers. Once a beat is detected, the function returns a non-zero delay indicating the QRS
	peak delay to the current sample.

 **********************************************************************************/

int16_t PT_StateMachine(int16_t datum)
{
	int16_t BeatDelay = 0;

	uint16_t PEAKI ;

	// ------- Preprocessing filtering and Peak detection --------- //
	LPFilter(&datum);										// LowPass filtering
	HPFilter();												// HighPass filtering

	PeakDtcBP(PT_dptr->HPF_val);							// Store BP signal highest peak
	
	DerivFilter();
	PeakDtcDR(PT_dptr->DRF_val);							// Store the highest slope for T wave discrimination

	SQRFilter();											//Squaring

	MVAFilter();
	PEAKI = PeakDtcI();

	// ---- Integrated Peak detection checks and blankTime ---- //
	if (!PEAKI && BlankTimeCnt)								// No beat, decrement BlankTime
	{
		if (--BlankTimeCnt == 0)							// If blanktime over place the oldest peak
			PEAKI = PEAKI_temp;
	}
	else if (PEAKI && !BlankTimeCnt)						// If no peak for peak for last 200msec, save the current peak
	{
		BlankTimeCnt = PT200MS;
		PEAKI_temp   = PEAKI;
		PEAKI = 0;
	}
	else if(PEAKI)											// If a bigger peak comes along, store it
	{
		if (PEAKI > PEAKI_temp)
		{
			BlankTimeCnt = PT200MS;
			PEAKI_temp = PEAKI;
			PEAKI = 0;
		}
		else if (--BlankTimeCnt == 0)
			PEAKI = PEAKI_temp;
		else
			PEAKI = 0;
	}

	// -- Run Different Phases of the Algo -> Learning Ph1, 2 and decision --//
	++Count_SinceRR;
	if (PT_dptr->PT_state == START_UP || PT_dptr->PT_state == LEARN_PH_1)		
	{ 
		if (PEAKI > 0)
			LearningPhase1(&PEAKI, &Best_PeakBP);
	}
	// ---- Once learning Phase 1 done, start storing beats ---- //
	else										
	{
		// ---- Is the peak taller than ThI1 and ThF1? ---- //
		if (PEAKI > PT_dptr->ThI1 && Best_PeakBP > PT_dptr->ThF1)
		{

			// ---- Initiated phase 2 ---- //
			if (PT_dptr->PT_state == LEARN_PH_2)
			{
				// ----- Update Integ & BP Th ------ //
				UpdateThI(&PEAKI, 0);
				UpdateThF(&Best_PeakBP, 0);

				// --- First RR interval --- //
				BeatDelay = GENERAL_DELAY + PT200MS;
				Count_SinceRR = 0;
				Old_PeakDR = Best_PeakDR;
				Best_PeakDR = 0;
				Best_PeakBP = 0;

				// --- Now we can compute RR intervals --- //
				PT_dptr->PT_state = DETECTING;

			}
			// ------ Learning phases are done! -------- //
			else
			{
			// --- T-Wave Test if RR < 360msec, is current slope lower 0.5prev_slope then noise --- //
				if (Count_SinceRR < PT360MS && (Best_PeakDR < (Old_PeakDR >> 2)))
				{
					// ----- Update Integ & BP Th ------ //
					UpdateThI(&PEAKI, 1);
					UpdateThF(&Best_PeakBP, 1);

				}
				else
				{
					// ----- Update Integ & BP Th && RR buffers ------ //
					UpdateThI(&PEAKI, 0);
					UpdateThF(&Best_PeakBP, 0);
					UpdateRR(Count_SinceRR);

					// --- Reset parameters --- //
					BeatDelay = GENERAL_DELAY + PT200MS;
					Count_SinceRR = 0;
					Old_PeakDR = Best_PeakDR;									// Store the derivative for T-wave test
					Best_PeakDR = Best_PeakBP = 0;

					SBcntI = 0;
					SB_peakBP = 0;
					SB_peakDR = 0;
					SB_peakI = 0;

				}
			}
		}
		// ------ If the peak is noise ------- //
		else if (PEAKI > 0)
		{
			// ----- Update Integ & BP Th ------ //
			UpdateThI(&PEAKI, 1);
			UpdateThF(&Best_PeakBP, 1);

			// ----- Store the peak for searchback ------ //
			if (PEAKI > SB_peakI && Count_SinceRR >= PT360MS)
			{
				SB_peakI = PEAKI;											// Store Integ Sig peak 
				SB_peakBP = Best_PeakBP;									// Store BP Sig peak
				SB_peakDR = Best_PeakDR;									// Derivative of SB point
				SBcntI = Count_SinceRR;										// Store Indice
			}

		}

	}

	// -- Do search-back if we have no beats in PT_dptr->RR_Missed_L -- //
	if (Count_SinceRR > PT_dptr->RR_Missed_L && SB_peakI > PT_dptr->ThI2 && PT_dptr->PT_state == DETECTING)
	{
		// ---- Checking the BP signal ---- //
		if (SB_peakBP > PT_dptr->ThF2)
		{
			// ----- Update Integ & BP Th && RR buffers ------ //
			UpdateThI(&SB_peakI, 0);
			UpdateThF(&SB_peakBP, 0);
			UpdateRR(SBcntI);

			// --- Reset parameters --- //
			BeatDelay = Count_SinceRR = Count_SinceRR - SBcntI;
			BeatDelay += (GENERAL_DELAY + PT200MS);
			Old_PeakDR = SB_peakDR;		// Store the derivative for T-wave test
			Best_PeakDR = Best_PeakBP = 0;

			SBcntI = 0;
			SB_peakBP = 0;
			SB_peakDR = 0;
			SB_peakI = 0;
		}
	}

	// ---- Emergency and Faulty Condition Reset ---- //
	// If algorithm doest not find a beat in 4sec, then it resets itself
	// and starts learning phases.
	if (Count_SinceRR > PT4000MS) {
		PT_init();
	}

	return (BeatDelay);
	
}


/**********************************************************************************

	Fuction Name: LearningPhase1

	Parameter:
	 Input	:	pkI - Pointer to the integrated signal peak
				pkBP- Pointer BP signal peak

	 Returns:	none - Updates the static values (st_mx_pk, st_mean_pk) and (st_mean_pkBP).

	Description: Computes the maximum peak in the past 2 seconds and also the mean of the
	peaks iteratively in both Integrated Signal and BP signal.

 **********************************************************************************/

void LearningPhase1(uint16_t *pkI, int16_t *pkBP)
{
	//---- Recursively compute the average and max of peaks ------ //
	if (*pkI > st_mx_pk) st_mx_pk = *pkI;

	// ---- If the very first time calling this function --- //
	if (PT_dptr->PT_state == START_UP) {
		PT_dptr->PT_state = LEARN_PH_1;
		st_mean_pk = *pkI;
		st_mean_pkBP = *pkBP; 
	}
	// ----- Continue averaging once still in learning ----- //
	else if(Count_SinceRR < PT2000MS){
		st_mean_pk = (st_mean_pk + *pkI) >> 1;
		st_mean_pkBP = (st_mean_pkBP + *pkBP) >> 1;
	}
	else {
		PT_dptr->PT_state = LEARN_PH_2;
		// ---- Integrated Signal Thresholds ------- //
		PT_dptr->SPKI = (st_mx_pk >> 1);
		PT_dptr->NPKI = (st_mean_pk >> 3);
		PT_dptr->ThI1 = PT_dptr->NPKI + ((PT_dptr->SPKI - PT_dptr->NPKI) >> 2);
		PT_dptr->ThI2 = PT_dptr->ThI1 >> 1;

		// -------- BP Signal Thresholds ---------- //
		PT_dptr->SPKF = (Best_PeakBP >> 1);
		PT_dptr->NPKF = (st_mean_pkBP >> 3);
		PT_dptr->ThF1 = PT_dptr->NPKF + ((PT_dptr->SPKF - PT_dptr->NPKF) >> 2);
		PT_dptr->ThF2 = PT_dptr->ThF1 >> 1;

	}
}

/**********************************************************************************

    Fuction Name: LPFilter

    Parameter:
     Input	:	none - Pointer to the input datum.

     Returns:	none - Updates the static value PT_dptr->LPF_val in place.

    Description: Low-pass filters the signal based on Pan-Tompkins Eq. 3,
	y[n] = 2*y[n-1] - y[n-2] + x[n] - 2 * x[n - 6] + x[n - 12] . This
	function implements the filter both in Direct Form I and II. Select the
	type employed by setting FILTER_FORM to 1 or 2. Delay of the filter is 5.

 **********************************************************************************/

void LPFilter(int16_t *val)
{
	// -- To avoid using modulo employ half-pointer -- //
	int16_t half_pointer, w;

	half_pointer = PT_dptr->LP_pointer - (LP_BUFFER_SIZE >> 1);

	if (half_pointer < 0) 
		half_pointer += LP_BUFFER_SIZE;

	

		// ------- Filter based on selected Form ------- //
#if (FILTER_FORM == 1)
		w = *val + (PT_dptr->LP_buf[1] << 1) - PT_dptr->LP_buf[0];
		*val = w - (PT_dptr->LP_buf[half_pointer] << 1) + PT_dptr->LP_buf[PT_dptr->LP_pointer];
		PT_dptr->LP_buf[PT_dptr->LP_pointer] = w;
#else
		w = (LP_y_old << 1) - LP_y_new + *val - (PT_dptr->LP_buf[half_pointer] << 1) + PT_dptr->LP_buf[PT_dptr->LP_pointer];
		LP_y_new = LP_y_old;
		LP_y_old = w;
		PT_dptr->LP_buf[PT_dptr->LP_pointer] = *val;
#endif
		// --- Avoid signal overflow by gaining down ---- //
		if (w >= 0)
			PT_dptr->LPF_val = w >> 5;
		else
			PT_dptr->LPF_val = (w >> 5) | 0xF800;

		if (++PT_dptr->LP_pointer == LP_BUFFER_SIZE) 
			PT_dptr->LP_pointer = 0;
}


/**********************************************************************************

Fuction Name: HPFilter

Parameter:
Input	:	none - Employs Current_Sample in place.

Returns:	none - Updates the static value PT_dptr->HPF_val  in place.

Description: High-pass filters the signal based on Pan-Tompkins Eq. 2.4 (Errata),
y[n] = y[n-1] + x[n-32]/32 - x[n]/32 + x[n - 16] - x[n - 17] . This
function implements the filter both in Direct Form I and II. Select the
type employed by setting FILTER_FORM to 1 or 2. Delay 16 samples.

**********************************************************************************/
void HPFilter(void)
{
	// -- To avoid using modulo employ half-pointer -- //
	int16_t half_pointer, h_prev_pointer;
	half_pointer = PT_dptr->HP_pointer - (HP_BUFFER_SIZE >> 1);

	if (half_pointer < 0)
		half_pointer += HP_BUFFER_SIZE;
	
	if (!half_pointer)
		h_prev_pointer = HP_BUFFER_SIZE - 1;
	else
		h_prev_pointer = half_pointer - 1;



	// ------- Filter based on selected Form ------- //
#if (FILTER_FORM == 1)
	y_h = PT_dptr->LPF_val + PT_dptr->HP_buf[0];
	PT_dptr->LPF_val = ((PT_dptr->HP_buf[PT_dptr->HP_pointer] - y_h) >> 5) + PT_dptr->HP_buf[half_pointer] - PT_dptr->HP_buf[h_prev_pointer];
	PT_dptr->HP_buf[PT_dptr->HP_pointer] = y_h;
#else
	y_h += (PT_dptr->HP_buf[PT_dptr->HP_pointer] >> 5) - (PT_dptr->LPF_val >> 5) + PT_dptr->HP_buf[half_pointer] - PT_dptr->HP_buf[h_prev_pointer];
	PT_dptr->HP_buf[PT_dptr->HP_pointer] = PT_dptr->LPF_val;
	
#endif
	// ------- Again slightly gaining down --------- //
	if (y_h >= 0)
		PT_dptr->HPF_val = (y_h >> 1);
	else
		PT_dptr->HPF_val = (y_h >> 1) | 0xF800;

	if (++PT_dptr->HP_pointer == HP_BUFFER_SIZE) PT_dptr->HP_pointer = 0;
}

/**********************************************************************************

Fuction Name: DerivFilter

Parameter:
Input	:	none - Employs Current_Sample.

Returns	:	none - Updates the value Current_Sample in place.

Description: Computes the signal derivative based on Pan-Tompkins Eq. 2.6 (Errata),
y[n] = 1/8(2x[n] + x[n - 1] - x[n - 3] - 2x[n - 4]) . Delay 2 samples.

**********************************************************************************/

void DerivFilter(void)
{
	// --- Since it is only a 5 point derivative filter we avoid using pointers and half pointers for further efficieny ---- //
	int16_t w;

	w = PT_dptr->DR_buf[0] - PT_dptr->DR_buf[2];
	w += ((PT_dptr->HPF_val - PT_dptr->DR_buf[3]) << 1);
	w >>= 3;
	PT_dptr->DR_buf[3] = PT_dptr->DR_buf[2];
	PT_dptr->DR_buf[2] = PT_dptr->DR_buf[1];
	PT_dptr->DR_buf[1] = PT_dptr->DR_buf[0];
	PT_dptr->DR_buf[0] = PT_dptr->HPF_val;
	PT_dptr->DRF_val = w;
}

/**********************************************************************************

Fuction Name: SQRFilter

Parameter:
Input	:	none - Employs Current_Sample.

Returns	:	none - Updates the value Current_Sample in place.

Description: Squares the signal based on Pan-Tompkins Eq. 10,
y[n] = x[n]^2. No delay.

**********************************************************************************/
void SQRFilter(void)
{
	// ------------ Avoiding Overflow -------------- //
	uint16_t temp;
	if (PT_dptr->DRF_val > SQR_LIM_VAL || PT_dptr->DRF_val < (-SQR_LIM_VAL))
		PT_dptr->SQF_val = UINT16_MAX;
	else
	{
		if (PT_dptr->DRF_val < 0)
			temp = (uint16_t)(-PT_dptr->DRF_val);
		else
			temp = (uint16_t)(PT_dptr->DRF_val);
		PT_dptr->SQF_val = temp*temp;
	}

	if (PT_dptr->SQF_val > SQR_LIM_OUT)
		PT_dptr->SQF_val = SQR_LIM_OUT;
}


/**********************************************************************************

Fuction Name: MVAFilter

Parameter:
Input	:	none - Employs Current_Sample.

Returns	:	none - Updates the value Current_Sample in place.

Description: Computes the rolling moving average of the input signal
based on Eq. 11 of Pan-Tompkins, y[n] = (1/N)[sum(x[1]+...+x[N])]. Delay 15 Samples.

**********************************************************************************/
void MVAFilter(void)
{
	//---- The MV_sum can easily overflow so we limit the bound by uint16 precision ------ //
	if (MV_sum < (UINT16_MAX - PT_dptr->SQF_val))
		MV_sum += PT_dptr->SQF_val;
	else
		MV_sum = UINT16_MAX;

	if (MV_sum > PT_dptr->MVA_buf[PT_dptr->MVA_pointer])
		MV_sum -= PT_dptr->MVA_buf[PT_dptr->MVA_pointer];
	else
		MV_sum = 0;

	PT_dptr->MVA_buf[PT_dptr->MVA_pointer] = PT_dptr->SQF_val;

	PT_dptr->MVA_val = MV_sum/(uint16_t) MVA_BUFFER_SIZE;

	if (PT_dptr->MVA_val > MVA_LIM_VAL)
		PT_dptr->MVA_val = MVA_LIM_VAL;

	if (++PT_dptr->MVA_pointer == MVA_BUFFER_SIZE) 
		PT_dptr->MVA_pointer = 0;
}


/**********************************************************************************

Fuction Name: PeakDtcI

Parameter:
Input	:	none - Employs Current_Sample.

Returns	:	p	 - The local maxima in intergrated signal.

Description: This is a simple peak detector for fiducial point detection in the integrated signal.
If the signal changes sign the value of the peak is asssumed ot be a peak.
if x[n-1] <= x[n] > x[n+1], then x[n] is a peak.

**********************************************************************************/
uint16_t PeakDtcI(void)
{
	uint16_t p;
	// ---------- Local maxima or not --------- //
	if (PT_dptr->MVA_val <= Prev_val && Prev_val > Prev_Prev_val) {
		p = Prev_val;
	}
	else {
		p = 0;
	}
	Prev_Prev_val = Prev_val;
	Prev_val = PT_dptr->MVA_val;

	return (p);
}

/**********************************************************************************

Fuction Name: PeakDtcDR

Parameter:
Input	:	none - Input from derivative filter

Returns	:	p	 - Stores the beat in Best_PeakDR for T-wave dsicrimation.

Description: This is a simple peak detector for fiducial point detection in Derivative Signal,
the strategy is to store the highest slope in the signal preceding the qrs so that if the next qrs is
suspicious to be a real beat could be compared against this. This helps us not to need to store a buffer.
For T-wave discrimination see T-wave identification of the paper.
if x[n-1] <= x[n] > x[n+1], then x[n] is a peak.

**********************************************************************************/
void PeakDtcDR(int16_t DR_sample)
{
	if (DR_sample < 0) DR_sample = -DR_sample;
	// ---------- Local maxima or not --------- //
	if (DR_sample <= Prev_valDR && Prev_valDR > Prev_Prev_valDR) {
		//-- For T-wave discrimination store the highest slope -- //
		if (Prev_valDR > Best_PeakDR) Best_PeakDR = Prev_valDR;
	}
	Prev_Prev_valDR = Prev_valDR;
	Prev_valDR = DR_sample;
}

/**********************************************************************************

Fuction Name: PeakDtcBP

Parameter:
Input	:	none - Input from BP signal.

Returns	:	none - Returns a local maxima in Best_PeakBP.

Description: This is a simple peak detector for fiducial point detection in BP Signal. 
Once a peak is detected in Integrated signal, the maximum peak in BP signal is also compared 
against adaptive thresholds.
if x[n-1] <= x[n] > x[n+1], then x[n] is a peak.

**********************************************************************************/
void PeakDtcBP(int16_t DR_sample)
{
	if (DR_sample < 0) DR_sample = -DR_sample;
	// ---------- Local maxima or not --------- //
	if (DR_sample <= Prev_valBP && Prev_valBP > Prev_Prev_valBP) {
		//-- For T-wave discrimination store the highest slope -- //
		if (Prev_valBP > Best_PeakBP) Best_PeakBP = Prev_valBP;
	}
	Prev_Prev_valBP = Prev_valBP;
	Prev_valBP = DR_sample;
}


/**********************************************************************************

Fuction Name: UpdateRR

Parameter:
Input	:	qrs		- The most recent RR interval.

Returns	:	none	- Updates the normal and selected RR values and updates the buffers.

Description: The most recent RR interval is passed to the function, if the RR interval is within
the normal range, RR_LOW_LIM, RR_HIGH_LIM, RR_MISSED_LIM, RR_AVRG2_buf is updated, otherwise
RR_AVRG1_buf is updated which is the average of the most recents beats. Note that two buffers
are employed for robust mean computation (Eq. 24-28). If the peak is irregular the thresholds
are reduced by 50%.

RR_Low_Lim		= 0.92*RR_M = ((92/100) * RR_M) = RR_M - (2/25)*RR_M
RR_High_Lim		= 1.16*RR_M = ((116/100) * RR_M) = RR_M + (4/25)*RR_M
RR_Missed_Lim	= 1.66*RR_M = ((166/100) * RR_M) = RR_M + (33/50)*RR_M

**********************************************************************************/
void UpdateRR(int16_t qrs)
{   
	// ---------- Update most 8 Recent RR mean Interval------------- //
	RR1_sum += qrs;
	RR1_sum -= PT_dptr->RR_AVRG1_buf[RR1_p];

	PT_dptr->RR_AVRG1_buf[RR1_p] = qrs;
	PT_dptr->Recent_RR_M = RR1_sum/RR_BUFFER_SIZE; 
	if (++RR1_p == RR_BUFFER_SIZE) 
		RR1_p = 0;



	// ------ Update Selected Beat RR mean if qrs in range --------- //
	if (qrs >= PT_dptr->RR_Low_L && qrs <= PT_dptr->RR_High_L) {
		// ------ Update selective RR mean ----- //
		RR2_sum += qrs;
		RR2_sum -= PT_dptr->RR_AVRG2_buf[RR2_p];

		PT_dptr->RR_AVRG2_buf[RR2_p] = qrs;
		PT_dptr->RR_M = RR2_sum / RR_BUFFER_SIZE;
		if (++RR2_p == RR_BUFFER_SIZE) 
			RR2_p = 0;

		// --------- Update Limits ------------ //
		PT_dptr->RR_Low_L = PT_dptr->Recent_RR_M - (PT_dptr->Recent_RR_M << 1) / 25;
		PT_dptr->RR_High_L = PT_dptr->Recent_RR_M + (PT_dptr->Recent_RR_M << 2) / 25;
		PT_dptr->RR_Missed_L = PT_dptr->RR_M + (PT_dptr->RR_M * 33) / 50;
		PT_dptr->HR_State = REGULAR_HR;
	}
	// -------- Irregular heart-rate ---------- //
	else {
		PT_dptr->RR_Missed_L = PT_dptr->Recent_RR_M + (PT_dptr->Recent_RR_M * 33) / 50;
		PT_dptr->ThI1 >>= 1;
		PT_dptr->ThF1 >>= 1;
		PT_dptr->HR_State = IRREGULAR_HR;
	}
	
}


/**********************************************************************************

Fuction Name: UpdateThI

Parameter:
Input	:	PEAKI	- Pointer to the overal peak
			NOISE_F	- flag that indicates whether the peak is noise or not.

Returns	:	none	- Recieves PEAKI and updates Signal peak estimate (SPKI), Noise peak
peak estimate (NPKI)

Description: This function recursively updates the adaptive noise and signal thresholds in
the Integrated signal. Implements Eq 12-16.

**********************************************************************************/
void UpdateThI( uint16_t *PEAKI, int8_t NOISE_F)
{
	// ------ Update Noise & Signal Estimate ------ //
	if (NOISE_F) {
		PT_dptr->NPKI -= PT_dptr->NPKI >> 3;
		PT_dptr->NPKI += *PEAKI >> 3;
	}
	else {
		PT_dptr->SPKI -= PT_dptr->SPKI >> 3;
		PT_dptr->SPKI += *PEAKI >> 3;
	}

	// --------- Update Thresholds ---------------- //
	PT_dptr->ThI1 = PT_dptr->NPKI + ((PT_dptr->SPKI - PT_dptr->NPKI) >> 2);
	PT_dptr->ThI2 = PT_dptr->ThI1 >> 1;
}


/**********************************************************************************

Fuction Name: UpdateThF

Parameter:
Input	:	PEAKF	- Pointer to the overal peak
			NOISE_F	- Flag that indicates whether the peak is noise or not.

Returns	:	none	- Recieves PEAKF and updates Signal peak estimate in BP signal(SPKF), Noise peak
peak estimate (NPKF)

Description: This function recursively updates the adaptive noise and signal thresholds in 
the BP signal. Implements Eq 17-21.

**********************************************************************************/
void UpdateThF(int16_t *PEAKF, int8_t NOISE_F)
{
	// ------ Update Noise & Signal Estimate ------ //
	if (NOISE_F) {
		PT_dptr->NPKF -= PT_dptr->NPKF >> 3;
		PT_dptr->NPKF += *PEAKF >> 3;
	}
	else {
		PT_dptr->SPKF -= PT_dptr->SPKF >> 3;
		PT_dptr->SPKF += *PEAKF >> 3;
	}

	// --------- Update Thresholds ---------------- //
	PT_dptr->ThF1 = PT_dptr->NPKF + ((PT_dptr->SPKF - PT_dptr->NPKF) >> 2);
	PT_dptr->ThF2 = PT_dptr->ThF1 >> 1;
}



/**************************************************
Helper functions for debugging and easy management.
One could use this to debug the algorithm in real-time or
query any parameter in the algorithm such as outputs of
the filters, thresholds, running heart-rate and etc.
***************************************************/

// ------Returns the state machine's state- SEE header for understanding underlying states ------ //
int16_t PT_get_State_output(void) {
	return (PT_dptr->PT_state);
}



// ------Returns LP filter value ------ //
int16_t PT_get_LPFilter_output(void) {
	return (PT_dptr->LPF_val);
}

// ------Returns HP filter value ------ //
int16_t PT_get_HPFilter_output(void) {
	return (PT_dptr->HPF_val);
}

// ------Returns Dr filter value ------ //
int16_t PT_get_DRFilter_output(void) {
	return (PT_dptr->DRF_val);
}

// ------Returns MVA filter value ------ //
uint16_t PT_get_MVFilter_output(void) {
	return (PT_dptr->MVA_val);
}

// ------Returns SQR filter value ------ //
uint16_t PT_get_SQRFilter_output(void) {
	return (PT_dptr->SQF_val);
}



/************************************
Returns instantanous heart rate per minute.

Input - Fs : Sampling Frequency of the signal
*************************************/
int16_t PT_get_ShortTimeHR_output(int16_t Fs) {
	return (60 / (PT_dptr->Recent_RR_M / Fs));
}

/************************************
Returns robust heart rate per minute

Input - Fs : Sampling Frequency of the signal
*************************************/
int16_t PT_get_LongTimeHR_output(int16_t Fs) {
	return (60 / (PT_dptr->RR_M / Fs));
}


// ------Returns the main threshold integrated signal Th value ------ //
uint16_t PT_get_ThI1_output(void) {
	return (PT_dptr->ThI1);
}

// ------Returns the main threshold BP signal Th value ------ //
int16_t PT_get_ThF1_output(void) {
	return (PT_dptr->ThF1);
}

// ------Returns Signal Level Estimate in Integrated Signal ----- //
uint16_t PT_get_SKPI_output(void) {
	return (PT_dptr->SPKI);
}

// ------Returns Noise Level Estimate in Integrated Signal ------ //
uint16_t PT_get_NPKI_output(void) {
	return (PT_dptr->NPKI);
}

// ------Returns Signal Level Estimate in BP Signal ------ //
int16_t PT_get_SPKF_output(void) {
	return (PT_dptr->SPKF);
}

// ------Returns Noise Level Estimate in BP Signal ------ //
int16_t PT_get_NPKF_output(void) {
	return (PT_dptr->NPKF);
}

// ------Returns HR state -> Regular:0, Irregular:1 ------ //
int16_t PT_get_HRState_output(void) {
	return (PT_dptr->HR_State);
}