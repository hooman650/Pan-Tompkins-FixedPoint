#ifndef _PANTOMPKINS_H_
#define _PANTOMPKINS_H_

#include <stdint.h>
#include <string.h>

/************************************************************
    Timing constants (assumption Sampling Frequency = 200 Hz)
 ************************************************************/
#define PT150MS				((int16_t)	(30))
#define PT200MS				((int16_t)	(40))
#define PT360MS				((int16_t)	(72))
#define PT1000MS			((int16_t)	(200))
#define PT2000MS			((int16_t)	(400))
#define PT4000MS			((int16_t)	(800))
#define GENERAL_DELAY		((int16_t)	(38))

/*************************************************************
	RR Limits constants for startup (92,116, 166 % OF 200)
*************************************************************/
#define RR92PERCENT			((int16_t)	(184))
#define RR116PERCENT		((int16_t)	(232))
#define RR166PERCENT		((int16_t)	(332))
/************************************************************
    General constants
 ************************************************************/
#define LP_BUFFER_SIZE               ((int16_t)  (12))
#define HP_BUFFER_SIZE               ((int16_t)  (32))
#define DR_BUFFER_SIZE               ((int16_t)  (4))
#define MVA_BUFFER_SIZE              ((int16_t)  (30))
#define RR_BUFFER_SIZE               ((int16_t)  (8))

/************************************************************
    PT constants
 ************************************************************/
#define FILTER_FORM		2
#define SQR_LIM_VAL			((int16_t)  (256))		// We have to limit the Squaring function to avoid overflow once squaring numbers.

#define SQR_LIM_OUT			((uint16_t) (30000))	// Hardlimiting output of Sqauring filter
#define MVA_LIM_VAL			((int16_t)	(32000))	// Limiting factor for MVA signal

 // values for State-Machine: PT_data.PT_state
#define START_UP		0
#define LEARN_PH_1		1
#define LEARN_PH_2		2
#define DETECTING		3



// Regular and Irregular Heart-rate
#define IRREGULAR_HR		1
#define REGULAR_HR			0



/************************************************************
    Data types
 ************************************************************/
struct PT_struct
{
	int16_t LP_pointer;
	int16_t HP_pointer;
	int16_t MVA_pointer;
	int16_t	PT_state;							//	State of the process
	int16_t Recent_RR_M;						//	Mean of most recent RR
	
	int16_t LPF_val;
	int16_t HPF_val;
	int16_t DRF_val;
	uint16_t SQF_val;
	uint16_t MVA_val;
	

	uint16_t ThI1;								// Threshold I1 (Integrated signal)
	uint16_t SPKI;								// Signal peak estimate Integrated
	uint16_t NPKI;								// Noise peak estimate Integrated
	uint16_t ThI2;								// Threshold I2 (Integrated signal)
	
	int16_t ThF1;								// Threshold F1 (Band-passed signal)
	int16_t SPKF;								// Signal peak estimate BP
	int16_t NPKF;								// Noise peak estimate BP
	int16_t ThF2;								// Threshold F2 (Band-passed signal)

	int16_t RR_M;								//  General mean of RR within acceptable range
	int16_t RR_Low_L;							//	RR Low limit 
	int16_t RR_High_L;							//	RR High limit 
	int16_t RR_Missed_L;						//	RR missed limit 
	int16_t HR_State;							//  HR-State can be regular or irregular

	int16_t LP_buf[LP_BUFFER_SIZE];				//  LP filter buffer
	int16_t HP_buf[HP_BUFFER_SIZE];				//  HP filter buffer
	int16_t DR_buf[DR_BUFFER_SIZE];				//  DR filter buffer
	uint16_t MVA_buf[MVA_BUFFER_SIZE];			//  MVA filter buffer
	int16_t RR_AVRG1_buf[RR_BUFFER_SIZE];		//  RR average 1 buffer
	int16_t RR_AVRG2_buf[RR_BUFFER_SIZE];		//  RR average 2 buffer
};

/**********************************************************************
    Function Prototypes
 **********************************************************************/
void PT_init(void);
int16_t PT_StateMachine(int16_t datum);
void LearningPhase1(uint16_t *pkI, int16_t *pkBP);
void LPFilter(int16_t *val);
void HPFilter(void);
void DerivFilter(void);
void SQRFilter(void);
void MVAFilter(void);
int16_t PeakDtcI(void);
void PeakDtcDR(int16_t DR_sample);
void PeakDtcBP(int16_t DR_sample);
void UpdateRR(int16_t qrs);
void UpdateThI(uint16_t *PEAKI, int8_t NOISE_F);
void UpdateThF(int16_t *PEAKF, int8_t NOISE_F);

/**********************************************************************
	Debuggin Functions
***********************************************************************/
int16_t PT_get_LPFilter_output(void);
int16_t PT_get_HPFilter_output(void);
int16_t PT_get_DRFilter_output(void);
uint16_t PT_get_MVFilter_output(void);
uint16_t PT_get_SQRFilter_output(void);
int16_t PT_get_ShortTimeHR_output(int16_t Fs);
int16_t PT_get_LongTimeHR_output(int16_t Fs);
uint16_t PT_get_ThI1_output(void);
int16_t PT_get_ThF1_output(void);
uint16_t PT_get_SKPI_output(void);
uint16_t PT_get_NPKI_output(void);
int16_t PT_get_SPKF_output(void);
int16_t PT_get_NPKF_output(void);
int16_t PT_get_HRState_output(void);

#endif

