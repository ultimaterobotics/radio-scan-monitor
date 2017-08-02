//#include <SDL/SDL.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

#include <sys/ipc.h>
#include <sys/shm.h>

#include <stdio.h>

#include "graph_tools.h"
#include "simplechart.h"
#include "csvReader.h"
#include "detector.h"

float wide_threshold = 2; //in dBm, difference between peak and background to start
float narrow_threshold = 10; //detection process

int debug_print = 0;

SDL_Texture* myTexture = NULL;

SDL_Surface *scr;
SDL_Window *screen;
SDL_Renderer* renderer = NULL;
SDL_Texture *scrt;
TTF_Font *font = NULL; 

void prepareOut(int w, int h)
{
//	if(debug_level >= 10) printf("prepareOut\n");
	SDL_Init(SDL_INIT_VIDEO);
	screen = SDL_CreateWindow("monitor",
                          SDL_WINDOWPOS_UNDEFINED,
                          SDL_WINDOWPOS_UNDEFINED,
                          w, h,
                          0);		
	renderer = SDL_CreateRenderer(screen, -1, 0);
	
	scrt = SDL_CreateTexture(renderer,
                               SDL_PIXELFORMAT_ARGB8888,
                               SDL_TEXTUREACCESS_STREAMING,
                               w, h);	
	SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);

	SDL_RenderClear(renderer);
	SDL_RenderPresent(renderer);	
}

void drawFrameM(uint8_t *drawBuf, int w, int h)
{
	SDL_UpdateTexture(scrt, NULL, drawBuf, w * 4);
	SDL_RenderClear(renderer);
	SDL_RenderCopy(renderer, scrt, NULL, NULL);
	SDL_RenderPresent(renderer);
}


float *full_spectrum_avg;
float *full_spectrum_max;
float *full_spectrum_avgZ;
float *full_frequencies;

float *full_detector_res;
float *full_detector_power;
float *full_detector_bw;
int detectors_count = 2;

int full_sp_size = 60000; //6GHz with 0.1MHz step
float full_sp_start_freq = 0; //in Hz
float full_sp_end_freq = 6000000000.0; //in Hz
float full_sp_freq_step;
int full_sp_min_filled_data = full_sp_size; //min > max indicating there is no data
int full_sp_max_filled_data = 0;

void init_spectrum()
{
	full_spectrum_avg = new float[full_sp_size];
	full_spectrum_avgZ = new float[full_sp_size];
	full_spectrum_max = new float[full_sp_size];
	full_frequencies = new float[full_sp_size];
	
	full_detector_res = new float[full_sp_size*detectors_count];
	full_detector_power = new float[full_sp_size*detectors_count];
	full_detector_bw = new float[full_sp_size*detectors_count];
	
	full_sp_freq_step = (full_sp_end_freq - full_sp_start_freq) / (float)full_sp_size;
	float cur_freq = full_sp_start_freq;
	for(int n = 0; n < full_sp_size; n++)
	{
		full_spectrum_avg[n] = -100; //default, never can be reached in hardware
		full_spectrum_avgZ[n] = 0.0000000001; //to avoid zero division in any case
		full_spectrum_max[n] = -100;
		full_frequencies[n] = cur_freq;
		cur_freq += full_sp_freq_step;
	}
	for(int n = 0; n < full_sp_size*detectors_count; n++)
	{
		full_detector_res[n] = 0;
		full_detector_power[n] = 0;
		full_detector_bw[n] = 0;
	}
}

int charts_size = 800;

float power_zoom = 100.0;

float rel_pos = 0.0;
float rel_length = 0.01;

CSimpleChart *main_chart;
CSimpleChart *main_chart_freq; //only for storing frequency data, not displaying them
CSimpleChart *zoom_chart;
CSimpleChart *zoom_chart_freq; //only for storing frequency data, not displaying them

CSimpleChart *detector_chart;
float zoom_f = 1;
float zero_level = -90;
int zoom_base_size = 2000;
int zoom_size = 2000;

void init_chart_detector()
{
	detector_chart = new CSimpleChart(800);
	detector_chart->setViewport(50, 10, 300, 200);
	detector_chart->setParameter("zero value", 0.0);
	detector_chart->setParameter("scale", 1.0);
}

void init_chart_main()
{
	if(main_chart != NULL)
	{
		delete main_chart;
		delete main_chart_freq;
	}
	main_chart_freq = new CSimpleChart(charts_size);
	main_chart = new CSimpleChart(charts_size);
	main_chart->setViewport(50, 10, 1000, 600);
	main_chart->setParameter("color", 255, 255, 255);
	main_chart->setParameter("draw axis", "yes");
	main_chart->setParameter("scaling", "manual");
	main_chart->setParameter("zero value", zero_level);
	main_chart->setParameter("scale", power_zoom);
}

void init_chart_zoom()
{
	if(zoom_chart != NULL)
	{
		delete zoom_chart;
		delete zoom_chart_freq;
	}
	zoom_chart_freq = new CSimpleChart(zoom_size);
	zoom_chart = new CSimpleChart(zoom_size);
	zoom_chart->setViewport(500, 10, 500, 300);
	zoom_chart->setParameter("color", 255, 255, 255);
	zoom_chart->setParameter("draw axis", "yes");
	zoom_chart->setParameter("scaling", "manual");
	zoom_chart->setParameter("zero value", zero_level);
	zoom_chart->setParameter("scale", power_zoom);
}

void rezoom()
{
	zoom_size = zoom_base_size * zoom_f;
	init_chart_zoom();
	zoom_chart->setParameter("scale", power_zoom);
	main_chart->setParameter("scale", power_zoom);
}

float mouse_rel_x = 0;

int fill_mode = 0; //0 - all data, 1 - only measured range

void fill_zoom_values()
{
	int start_idx = 0;
	
	if(fill_mode == 1)
		start_idx = full_sp_min_filled_data;
	
	zoom_chart->clear();
	int rbg = start_idx + (charts_size - zoom_size) * mouse_rel_x;
	int red = rbg + zoom_size;
	if(red >= full_sp_size) red = full_sp_size;
	for(int r = rbg; r < red; r ++)
	{
		if(full_spectrum_avgZ[r] > 0.5)
//			zoom_chart->addV(full_spectrum_avg[r] / full_spectrum_avgZ[r]);
			zoom_chart->addV(full_spectrum_max[r]);
		else 
			zoom_chart->addV(-100);
		zoom_chart_freq->addV(full_frequencies[r]);			
	}
}

void fill_spectrum_data()
{
	if(fill_mode == 0)
		charts_size = full_sp_size;
	else 
		charts_size = full_sp_max_filled_data - full_sp_min_filled_data;
	if(charts_size <= 0) return;
	init_chart_zoom();
	init_chart_main();
	main_chart->clear();
	int start_idx = 0;
	int end_idx = full_sp_size;
	
	if(fill_mode == 1)
	{
		start_idx = full_sp_min_filled_data;
		end_idx = full_sp_max_filled_data;
	}
	
	for(int r = start_idx; r < end_idx; r ++)
	{
		if(full_spectrum_avgZ[r] > 0.5)
//			main_chart->addV(full_spectrum_avg[r] / full_spectrum_avgZ[r]);
			main_chart->addV(full_spectrum_max[r]);
		else 
			main_chart->addV(-100);
		main_chart_freq->addV(full_frequencies[r]);			
	}
	fill_zoom_values();
}

#define SHM_SIZE 1000000
//1mb - more than we normally need, just in case of strange input parameters

uint8_t *shared_memory = NULL;

void shared_mem_init()
{
    key_t key = 47192032; //some random number that must be the same in gr-scan shared mem module
    int shmid;

    if ((shmid = shmget(key, SHM_SIZE, IPC_CREAT | 0666)) < 0) {
        printf("shmget error!\n");
        exit(1);
    }

    if ((shared_memory = (uint8_t*)shmat(shmid, NULL, 0)) == (uint8_t *) -1) {
        printf("shmat error!\n");
        exit(1);
    }
}

int need_update_detector = 0;

void load_scan_data()
{
	float *f_shm = (float*)shared_memory;
	int *i_shm = (int*)shared_memory;
//	float start_freq = f_shm[1];
//	float end_freq = f_shm[2];
	float gain_mod = 0;//f_shm[3];
	float gain_add = 10.0;
	
	float norm_avg_param = 0.97;
	float max_mult_param = 1.01;
	int num_points = i_shm[4];
	
//	float len = num_points;
	
	float centr_freq = f_shm[5 + num_points];
	int centr_pos = (centr_freq - full_sp_start_freq) / full_sp_freq_step;
	
	int fill_cp = (full_sp_min_filled_data + full_sp_max_filled_data)/2;
	if(centr_pos - fill_cp < 2000 && !need_update_detector) need_update_detector = 1;
	
	for(int r = num_points/4; r < 3*num_points/4; r++)
	{
		float freq = f_shm[5 + r*2];
		float value = f_shm[6 + r*2];
		if(gain_mod != 0) value += gain_add;
		if(r < 10 || r > num_points-10 || (r > num_points/2-4 && r < num_points/2+4)) continue;

		int freq_pos = (freq - full_sp_start_freq) / full_sp_freq_step;
		if(freq_pos < 1 || freq_pos >= full_sp_size) continue;
		if(freq_pos < full_sp_min_filled_data) full_sp_min_filled_data = freq_pos;
		if(freq_pos > full_sp_max_filled_data) full_sp_max_filled_data = freq_pos;
		if(full_spectrum_avgZ[freq_pos] < 1)
		{
			full_spectrum_avg[freq_pos] = value;
			full_spectrum_avgZ[freq_pos] = 1.0;
			full_spectrum_max[freq_pos] = value;
		}
		else
		{
			full_spectrum_avg[freq_pos] *= norm_avg_param;
			full_spectrum_avg[freq_pos] += (1.0 - norm_avg_param) * value;
			full_spectrum_avgZ[freq_pos] = 1.0;
			full_spectrum_max[freq_pos] *= max_mult_param;
			if(value > full_spectrum_max[freq_pos])
				full_spectrum_max[freq_pos] = value;
		}
	}
	fill_spectrum_data();
}

int last_scan_id = 0; //last scan ID that we processed

void shared_mem_controller()
{
	if(!shared_memory)
	{
		shared_mem_init();
		last_scan_id = ((int*)shared_memory)[0];
	}
	int *i_shm = (int*)shared_memory;
	if(i_shm[0] != last_scan_id)
	{
		load_scan_data();
		last_scan_id = i_shm[0];
	}
}

void clear_all_data()
{
	int *i_shm = (int*)shared_memory;
	i_shm[0] = 0;
	for(int x = 0; x < full_sp_size; x++)
	{
		full_spectrum_avg[x] = -100;
		full_spectrum_avgZ[x] = 0.0000001;
		full_spectrum_max[x] = -100;
	}
	full_sp_max_filled_data = 0;
	full_sp_min_filled_data = full_sp_size;
}


void draw_charts(uint8_t *draw_pix, int w, int h)
{
	main_chart->draw(draw_pix, w, h);
	zoom_chart->draw(draw_pix, w, h);

	int zbx = zoom_chart->getX();
	int zex = zoom_chart->getX() + zoom_chart->getSizeX();
	int zby = zoom_chart->getY();
	int zey = zoom_chart->getY() + zoom_chart->getSizeY();

	int mbx = main_chart->getX();
	int mex = main_chart->getX() + main_chart->getSizeX();
	int mby = main_chart->getY();
	int mey = main_chart->getY() + main_chart->getSizeY();
	
	//drawing axes:
	//drawing grids:
	for(float ygrid = -90; ygrid <= 0; ygrid += 10)
	{
		int gy = main_chart->getValueY(ygrid);
		if(gy < mby || gy >= mey) continue;
		if(gy < 0 || gy >= h) continue;
		if(mbx < 0) mbx = 0;
		if(mex > w) mex = w;
		for(int x = mbx; x < mex; x++)
		{
			if(x >= zbx && x <= zex && gy >= zby && gy <= zey) continue;
			((unsigned int*)draw_pix)[gy*w + x] = 0xFFFFFF;
		}
	}
	for(float ygrid = -90; ygrid <= 0; ygrid += 10)
	{
		int gy = zoom_chart->getValueY(ygrid);
		if(gy < zby || gy >= zey) continue;
		if(gy < 0 || gy >= h) continue;
		if(zbx < 0) zbx = 0;
		if(zex > w) zex = w;
		for(int x = zbx; x < zex; x++)
			((unsigned int*)draw_pix)[gy*w + x] = 0xFFFFFF;
	}
	
//	detector_chart->draw(draw_pix, w, h);
}


uint8_t bmp_s[2];
uint32_t bmp_sign[3];
uint32_t bmp_header[11];
typedef struct sDIBHeader {
	uint32_t head_size; //12
	uint16_t width;
	uint16_t height;
	uint16_t planes;
	uint16_t bitpp;
}sDIBHeader;
sDIBHeader dbhead;

void init_bmp_header()
{	
	bmp_s[0] = 'B';
	bmp_s[1] = 'M';
	bmp_sign[0] = 0; //size
	bmp_sign[1] = 0;
	bmp_sign[2] = 14+12; //offset
	
}


void store_bmp_file(int w, int h, uint8_t *buf)
{
	char fname[1024];
	timeval curTime;
	gettimeofday(&curTime, NULL);

	init_bmp_header();
	sprintf(fname, "scan_res_%ld.bmp", curTime.tv_sec);
	int fl = open(fname, O_WRONLY | O_CREAT, 0b110110110);
	int wr1 = write(fl, bmp_s, 2);
	if(wr1 < 1) if(debug_print) printf("bmp write error\n");
	bmp_sign[0] = w*h*4 + 14+12;
	wr1 = write(fl, bmp_sign, 12);
	if(wr1 < 1) if(debug_print) printf("bmp write error\n");
	dbhead.head_size = 12;
	dbhead.width = w;
	dbhead.height = h;
	dbhead.planes = 1;
	dbhead.bitpp = 32;
	wr1 = write(fl, &dbhead, sizeof(dbhead));
	if(wr1 < 1) if(debug_print) printf("bmp write error\n");
	for(int ln = 0; ln < h; ln++)
	{
		wr1 = write(fl, buf+(h-1-ln)*w*4, w*4);
		if(wr1 < 1) if(debug_print) printf("bmp write error\n");
	}
	close(fl);
}

void save_logs()
{
	char logfn_full[512];
	char logfn_short[512];
	
	time_t rawtime;
	time (&rawtime);
	struct tm * curTm = localtime(&rawtime);
	sprintf(logfn_full, "full_scan_y%d_m%d_d%d_h%d_m%d_s%d.csv", (2000+curTm->tm_year-100), curTm->tm_mon, curTm->tm_mday, curTm->tm_hour, curTm->tm_min, curTm->tm_sec);
	sprintf(logfn_short, "short_scan_y%d_m%d_d%d_h%d_m%d_s%d.csv", (2000+curTm->tm_year-100), curTm->tm_mon, curTm->tm_mday, curTm->tm_hour, curTm->tm_min, curTm->tm_sec);
	
	int log_full = open(logfn_full, O_WRONLY | O_CREAT, 0b110110110);
	char line[1024];
	
	for(int r = 1; r < full_sp_size; r++)
	{
		int lng = 0;
		if(full_spectrum_avgZ[r] > 0.5)
			lng = sprintf(line, "%g;%g;%g\n", full_frequencies[r], full_spectrum_avg[r] / full_spectrum_avgZ[r], full_spectrum_max[r]);
		else
			lng = sprintf(line, "%g;%g;%g\n", full_frequencies[r], -100.0, -100.0);
		int wlen = write(log_full, line, lng);
		if(wlen != lng)
			if(debug_print) printf("write error\n");
	}
	close(log_full);		

	int log_short = open(logfn_short, O_WRONLY | O_CREAT, 0b110110110);
	for(int r = full_sp_min_filled_data+1; r < full_sp_max_filled_data; r++)
	{
		int lng = 0;
		if(full_spectrum_avgZ[r] > 0.5)
			lng = sprintf(line, "%g;%g;%g\n", full_frequencies[r], full_spectrum_avg[r] / full_spectrum_avgZ[r], full_spectrum_max[r]);
		else continue;
		int wlen = write(log_short, line, lng);
		if(wlen != lng)
			if(debug_print) printf("write error\n");
	}
	close(log_short);		
}

sSignalDetector *detectors;

void init_detectors()
{
	detectors = new sSignalDetector[detectors_count];
	sprintf(detectors[0].name, "wide");
	detectors[0].single_peak = 1;
	detectors[0].standard_max_frequency_MHz = 99999;
	detectors[0].standard_min_frequency_MHz = 0;
	detectors[0].center_width_kHz = 7000;
	detectors[0].left_shift_kHz = 7000;
	detectors[0].left_width_kHz = 2000;
	detectors[0].right_shift_kHz = 7000;
	detectors[0].right_width_kHz = 2000;
	detectors[0].left_relative_power = 0.0;
	detectors[0].right_relative_power = 0.0;

	sprintf(detectors[1].name, "narrow");
	detectors[1].single_peak = 1;
	detectors[1].standard_max_frequency_MHz = 99999;
	detectors[1].standard_min_frequency_MHz = 0;
	detectors[1].center_width_kHz = 300;
	detectors[1].left_shift_kHz = 1550;
	detectors[1].left_width_kHz = 1400;
	detectors[1].right_shift_kHz = 1550;
	detectors[1].right_width_kHz = 1400;
	detectors[1].left_relative_power = 0.9;
	detectors[1].right_relative_power = 0.9;
}

float freq_step_hz = 100000;

typedef struct sDetectedSignal
{
	int type;
	float central_frequency;
	float BW;
	float power;
}sDetectedSignal;
#define MAX_DETECTIONS 10000
sDetectedSignal detected_signals[MAX_DETECTIONS];
int detected_signals_count;

void clear_detected_signals()
{
	detected_signals_count = 0;
}
void add_detected_signal(int type, float center, float BW, float power)
{
	if(detected_signals_count >= MAX_DETECTIONS) return;
	for(int s = 0; s < detected_signals_count; s++)
	{
		if(detected_signals[s].type != type) continue;
		float c2 = detected_signals[s].central_frequency;
		float bw2 = detected_signals[s].BW;
		if((center - BW < c2 && center + BW > c2) || (c2 - bw2 < center && c2 + bw2 > center)) //intersect
		{
			if(power > detected_signals[s].power)
			{
				detected_signals[s].central_frequency = center;
				detected_signals[s].BW = BW;
				detected_signals[s].power = power;
				return;
			}
			else
				return;
		}
	}
	detected_signals[detected_signals_count].type = type;
	detected_signals[detected_signals_count].central_frequency = center;
	detected_signals[detected_signals_count].BW = BW;
	detected_signals[detected_signals_count].power = power;
	detected_signals_count++;
}
void print_detected_signals()
{
	for(int s = 0; s < detected_signals_count; s++)
//	if((detected_signals[s].type == 0 && detected_signals[s].BW > 4.5 ) || (detected_signals[s].type == 1 && detected_signals[s].BW < 5) )
	{
		printf("\ttype: %s F %.1f MHz P %.0f dBm BW %.1f MHz\n", detectors[detected_signals[s].type].name, detected_signals[s].central_frequency, detected_signals[s].power, detected_signals[s].BW);
	}
}

void run_detectors()
{
	clear_detected_signals();
	int has_detected = 0;
	detector_chart->clear();
	for(int d = 0; d < detectors_count; d++)
	{
		int dwidth = detectors[d].get_window_width_points(freq_step_hz);
//		float df_min = detectors[d].get_min_freq();
//		float df_max = detectors[d].get_max_freq();
		for(int x = full_sp_min_filled_data; x < full_sp_max_filled_data - dwidth; x++)
		{
			float res_power = 0;
			float res_bw = 0;
			float det_level = detectors[d].process_data(full_spectrum_max + x, full_frequencies[x + dwidth/2], freq_step_hz, dwidth, &res_power, &res_bw);
			int sp_idx = x + dwidth/2;
			full_detector_res[full_sp_size*d + sp_idx] = det_level;
			full_detector_power[full_sp_size*d + sp_idx] = res_power;
			full_detector_bw[full_sp_size*d + sp_idx] = res_bw;
//			printf("%g (%g - %g);%g;%g;%g\n", full_frequencies[x + dwidth/2], full_frequencies[x], full_frequencies[x+dwidth], det_level, res_power, res_bw);
			if(0)if(fabs(full_frequencies[x + dwidth/2] - 2412500000) < 110000)
			{
				printf("Freq: %g, dwidth %d, dlvl %g, dpwr %g\n", full_frequencies[x + dwidth/2], dwidth, det_level, res_power);
				for(int n = x + dwidth/2 - 100; n < x + dwidth/2 + 100; n++)
					printf("%g\n", full_spectrum_max[n]);
			}
			if(full_frequencies[x + dwidth/2] > 1200*1000000.0 && full_frequencies[x + dwidth/2] < 1265*1000000.0)
			{
//				printf("%g;%g;%g;%g\n", full_frequencies[x + dwidth/2], det_level, res_power, res_bw);
				if(d==1)detector_chart->addV(det_level);
			}
		}
	}
	int loc_max_pos = 0;
	int loc_max_length_right = 0;
	int loc_max_length_left = 0;
	int min_peak_length = 3; //don't even consider spikes as signal
	int valid_peak_length = 6;
	int valid_max_length = 300; //30MHz, anything more won't be considered as signal
	float min_threshold = wide_threshold; //anything less won't be considered as a signal at all
	float in_peak_threshold = 0.7;
	float loc_max = 0;
//	float loc_max_freq = 0;
//	float loc_start_bw = 0;
//	float loc_end_bw = 0;
	for(int d = 0; d < detectors_count; d++)
	{
		if(d == 1) min_threshold = narrow_threshold; //strict rules for narrow filters
		for(int x = full_sp_min_filled_data; x < full_sp_max_filled_data; x++)
	{
		
		int didx = full_sp_size*d;
		float dlvl = full_detector_res[didx + x];
		if(dlvl > loc_max && full_detector_res[didx + x] > min_threshold)
		{
			loc_max = full_detector_res[didx + x];
			loc_max_pos = x;
//				loc_max_freq = full_frequencies[x];
			loc_max_length_right = 0;
			loc_max_length_left = 0;
		}
		if(dlvl > in_peak_threshold * loc_max)
		{
			loc_max_length_right++;
		}
		else
		{
			if(loc_max_length_right > min_peak_length)
			{
				for(int y = loc_max_pos; y > 0; y--)
				{
					if(full_detector_res[didx + y] < in_peak_threshold * loc_max) break;
					loc_max_length_left++;
				}
				int lng = loc_max_length_left + loc_max_length_right;
				if(lng >= valid_peak_length && lng < valid_max_length)
				{
					float center_freq = 0.5*(full_frequencies[loc_max_pos + loc_max_length_right] + full_frequencies[loc_max_pos - loc_max_length_left]) * 0.000001;
					float sig_power = full_detector_power[didx + loc_max_pos];
					float sig_bw = (full_frequencies[loc_max_pos + loc_max_length_right] - full_frequencies[loc_max_pos - loc_max_length_left])* 0.000001;
					float avg_power = 0;
					float max_power = -100;
					float min_power = 0;
					float avgZ = 0;
					for(int px = loc_max_pos - loc_max_length_left; px < loc_max_pos + loc_max_length_right; px++)
					{
						if(full_spectrum_max[px] > max_power) max_power = full_spectrum_max[px];
						if(full_spectrum_max[px] < min_power) min_power = full_spectrum_max[px];
						avg_power += full_spectrum_max[px];
						avgZ++;
					}
					avg_power /= avgZ;

					sig_power = 0.7*avg_power + 0.3*max_power;

					int bw_start = loc_max_pos - loc_max_length_left, bw_end = loc_max_pos + loc_max_length_right;
					float sp_avgs = avg_power * 0.9;
					for(int px = loc_max_pos; px < loc_max_pos + loc_max_length_right; px++)
					{
						sp_avgs *= 0.3;
						sp_avgs += 0.7*full_spectrum_max[px];
						if(sp_avgs < avg_power*1.05)
						{
							bw_end = px;
							break;
						}
					}
					sp_avgs = avg_power * 1.1;
					for(int px = loc_max_pos+1; px > loc_max_pos - loc_max_length_left; px--)
					{
						sp_avgs *= 0.3;
						sp_avgs += 0.7*full_spectrum_max[px];
						if(sp_avgs < avg_power*1.05)
						{
							bw_start = px;
							break;
						}
					}
					
					avg_power = 0;
					max_power = -100;
					min_power = 0;
					avgZ = 0;
					for(int px = bw_start; px < bw_end; px++)
					{
						if(full_spectrum_max[px] > max_power) max_power = full_spectrum_max[px];
						if(full_spectrum_max[px] < min_power) min_power = full_spectrum_max[px];
						avg_power += full_spectrum_max[px];
						avgZ++;
					}
					avg_power /= avgZ;
					sig_power = 0.7*avg_power + 0.3*max_power;
					
					float below_avg = 0;
					for(int px = bw_start; px < bw_end; px++)
						if(full_spectrum_max[px] < avg_power) below_avg += 1.0;
					below_avg /= (float)(bw_end - bw_start);
					
					sig_bw = (full_frequencies[bw_end] - full_frequencies[bw_start])*0.000001;
					center_freq = 0.5*(full_frequencies[bw_end] + full_frequencies[bw_start]) * 0.000001;
					if(sig_bw > 0.5) sig_bw -= 0.2; //compensate for edge effects for too narrow bands
					if(sig_bw > 0.6) sig_bw -= 0.2;

					double pow_integr = 0.00000001;
					double pow_high_05 = 0.00000001;
					for(int px = bw_start; px < bw_end; px++)
					{
						pow_integr += pow(10.0, full_spectrum_max[px] * 0.1);
						if(full_spectrum_max[px] > min_power + 0.5*(max_power - min_power))
							pow_high_05 += pow(10.0, full_spectrum_max[px] * 0.1);
					}

					sig_power = 10.0 * log10(pow_integr);
					float sig_power_high_05 = 10.0 * log10(pow_high_05);
					
					if(sig_power_high_05 < sig_power-3) continue; //too unequal power

// 						if(sig_power > -80.0)
					if(!has_detected)
					{
						printf("\nDetected signals:\n");
						has_detected = 1;
					}
					//printf("F %.1f P %.0f BW %.1f sdv %g mm %g\n", center_freq, sig_power, sig_bw, power_sdv, max_power - min_power);
					if(sig_power > -70 && below_avg < 0.7)
					if((d == 0 && sig_bw > 4.5 ) || (d == 1 && sig_bw < 5) )
					{
						add_detected_signal(d, center_freq, sig_bw, sig_power);
						//printf("\ttype: %s F %.1f MHz pow %g pow09 %g\n", detectors[d].name, center_freq, sig_power, sig_power_high_05);
						
					}
//							printf("\ttype: %s F %.1f MHz P %.0f dBm BW %.1f MHz\n", detectors[d].name, center_freq, sig_power, sig_bw);
				}
			}
			loc_max = 0;
			loc_max_length_left = 0;
			loc_max_length_right = 0;
		}
	}
	
	}
	if(!has_detected)
		printf("\n no signals detected \n");
	else
		print_detected_signals();
}

int main(int argc, char* argv[])
{
	if(debug_print) printf("starting:\n");
	init_spectrum();
	if(debug_print) printf("memory allocated\n");
	
	SDL_Surface *msg = NULL;
	
	int w = 1100;//1000;
	int h = 680;//500;
	
	uint8_t *drawPix = (uint8_t*)malloc(w*h*4);
	prepareOut(w, h);
	if(debug_print) printf("SDL created\n");
	//text output
	SDL_Color textColor = { 255, 255, 255 }; 	
	TTF_Init();
	font = TTF_OpenFont( "../Ubuntu-R.ttf", 14 ); //filename, wanted font size
	if(font == NULL ) 
	{
		font = TTF_OpenFont( "Ubuntu-R.ttf", 14 ); //filename, wanted font size
		if(font == NULL ) 
			printf("can't load font!\n");
	}
	if(debug_print) printf("font ok\n");
	
	init_detectors();
	init_chart_detector();
	
	int done = 0;
	
	timeval curTime, prevTime, zeroTime, detReqTime;
	gettimeofday(&prevTime, NULL);
	gettimeofday(&zeroTime, NULL);
	
	int mbtn_l, mbtn_r, mbtn_m;
	int shift_pressed = 0, alt_pressed = 0;
	if(debug_print) printf("starting main cycle\n");
	while( !done ) 
	{ 

//===========INPUT PROCESSING==================================
		int mouse_x = 0, mouse_y = 0;
		SDL_Event event;
		int store_image = 0;
		if(debug_print) printf("polling events:\n");
		while( SDL_PollEvent( &event ) ) 
		{ 
			if( event.type == SDL_QUIT ) 
			{
				done = 1; 
				printf("SDL quit message\n");
			} 
			if( event.type == SDL_KEYDOWN ) 
			{
				if(event.key.keysym.scancode == SDL_SCANCODE_ESCAPE) 
				{ 
					done = 1;
					printf("esc pressed\n");
				} 
				if(event.key.keysym.scancode == SDL_SCANCODE_SPACE) 
				{ 
					store_image = 1;
				} 
				if(event.key.keysym.scancode == SDL_SCANCODE_RETURN || event.key.keysym.scancode == SDL_SCANCODE_RETURN2 || event.key.keysym.scancode == SDL_SCANCODE_KP_ENTER) 
				{
					save_logs();
				}
				if(event.key.keysym.scancode == SDL_SCANCODE_F) 
				{
					fill_mode = !fill_mode;
					fill_spectrum_data();
//					printf("%d  -  %d, fill mode %d\n", full_sp_min_filled_data, full_sp_max_filled_data, fill_mode);
				}
				if(event.key.keysym.scancode == SDL_SCANCODE_R) 
				{
					clear_all_data();
					fill_spectrum_data();
				}
				if(event.key.keysym.scancode == SDL_SCANCODE_D) 
				{
					run_detectors();
				}
				if(event.key.keysym.scancode == SDL_SCANCODE_UP) 
				{ 
					if(shift_pressed)
						power_zoom *= 0.9;
					else if(alt_pressed)
						zero_level--;
					else 
						zoom_f *= 0.8;
					rezoom();
				} 
				if(event.key.keysym.scancode == SDL_SCANCODE_DOWN) 
				{ 
					if(shift_pressed)
						power_zoom *= 1.1;
					else if(alt_pressed)
						zero_level++;
					else 
						zoom_f *= 1.2;
					rezoom();
				} 

			}
			if(event.type == SDL_MOUSEWHEEL)
			{
				if(event.wheel.y > 0)
				{
					if(shift_pressed)
						power_zoom *= 0.9;
					else if(alt_pressed)
						zero_level--;
					else 
						zoom_f *= 0.8;
					rezoom();
				}
				if(event.wheel.y < 0)
				{
					if(shift_pressed)
						power_zoom *= 1.1;
					else if(alt_pressed)
						zero_level++;
					else 
						zoom_f *= 1.2;
					rezoom();
				}
			}
		} 
		if(debug_print) printf("get mouse state:\n");

//==============MOUSE==========================================
		uint8_t mstate;
		//mstate = SDL_GetRelativeMouseState (&mx, &my);
		mstate = SDL_GetMouseState(&mouse_x, &mouse_y);

		int lclk = 0;
		if (mstate & SDL_BUTTON (1))
		{
			if(mbtn_l == 0) lclk = 1;
			mbtn_l = 1;
		}
		else
		{
			mbtn_l = 0;
		}

		if (mstate & SDL_BUTTON (2)) mbtn_m = 1;
		if (mstate & SDL_BUTTON (3)) mbtn_r = 1;
		
		mbtn_m = mbtn_r + lclk;
		mbtn_r = mbtn_m;
		if(mouse_y > h*0.7)
		{
			if(main_chart != NULL)
				if(main_chart->getSizeX() > 1)
				{
					mouse_rel_x = (float)(mouse_x - main_chart->getX()) / (float)main_chart->getSizeX();
					if(mouse_rel_x < 0) mouse_rel_x = 0;
					if(mouse_rel_x > 1) mouse_rel_x = 1;
				}
		}
		
//==============MOUSE END======================================

//==============KEYBOARD==========================================
		if(debug_print) printf("get keyboard state:\n");
		uint8_t *keys = (uint8_t *)SDL_GetKeyboardState( NULL );
		
		//All key codes:
		//http://sdljava.sourceforge.net/docs/api/sdljava/x/swig/SDLKeyValues.html
		
		if(keys[SDL_SCANCODE_ESCAPE])
			done = 1;
		shift_pressed = keys[SDL_SCANCODE_LSHIFT];
		alt_pressed = keys[SDL_SCANCODE_LALT];

		if(keys[SDL_SCANCODE_Q])
			detectors[1].center_width_kHz *= 1.001;
		if(keys[SDL_SCANCODE_A])
			detectors[1].center_width_kHz *= 0.999;
		if(keys[SDL_SCANCODE_W])
		{
			detectors[1].left_shift_kHz *= 1.001;
			detectors[1].right_shift_kHz *= 1.001;
		}
		if(keys[SDL_SCANCODE_S])
		{
			detectors[1].left_shift_kHz *= 0.999;
			detectors[1].right_shift_kHz *= 0.999;
		}
		if(keys[SDL_SCANCODE_Z])
		{
			detectors[1].left_width_kHz *= 1.001;
			detectors[1].right_width_kHz *= 1.001;
		}
		if(keys[SDL_SCANCODE_X])
		{
			detectors[1].left_width_kHz *= 0.999;
			detectors[1].right_width_kHz *= 0.999;
		}
		

//==============KEYBOARD END======================================
		
//===========INPUT PROCESSING END==============================

		if(debug_print) printf("input processed\n");

		if(need_update_detector == 1)
		{
			gettimeofday(&detReqTime, NULL);
			need_update_detector = 2;
		}
		gettimeofday(&curTime, NULL);
		int dT = (curTime.tv_sec - prevTime.tv_sec) * 1000000 + (curTime.tv_usec - prevTime.tv_usec);

		if(need_update_detector == 2)
		{
			int det_dt = (curTime.tv_sec - detReqTime.tv_sec) * 1000000 + (curTime.tv_usec - detReqTime.tv_usec);
			if(det_dt > 1000000)
			{
				need_update_detector = 0;
				run_detectors();
			}
		}

		memset(drawPix, 0, w*h*4);
		prevTime = curTime;
		
//		rel_pos += 0.0005;
		if(debug_print) printf("shared memory controller:\n");
		shared_mem_controller();
		if(debug_print) printf("done\n");
		if(full_sp_max_filled_data <= full_sp_min_filled_data)
		{ 
			if(debug_print) printf("no data, ending cycle\n");
			SDL_RenderClear(renderer);
			SDL_RenderPresent(renderer);
			continue;
		}
		if(debug_print) printf("filling spectrum\n");
		fill_spectrum_data();
//		fill_zoom_values();
		if(debug_print) printf("drawing charts\n");
		draw_charts(drawPix, w, h);
		
		if(store_image)
		{
			if(debug_print) printf("storing image\n");
			store_image = 0;
			store_bmp_file(w, h, drawPix);
		}
		
		if(debug_print) printf("rendering\n");
		SDL_RenderClear(renderer);
		SDL_UpdateTexture(scrt, NULL, drawPix, w * 4);
		SDL_RenderCopy(renderer, scrt, NULL, NULL);

		int curY = 25, curDY = 15;
		char outstr[128];
		SDL_Rect mpos;
		SDL_Point cpos;

		static float fps = 0;

		if(debug_print) printf("text out start\n");
		
		float mDT = 1000000.0 / (float)dT;
		fps = fps*0.9 + 0.1*mDT;
		sprintf(outstr, "fps %.0f cw %g lf %g lw %g", fps, detectors[1].center_width_kHz, detectors[1].left_shift_kHz, detectors[1].left_width_kHz);
		msg = TTF_RenderText_Solid(font, outstr, textColor); 
		mpos.x = 5; mpos.y = curY; curY += curDY;
		mpos.w = msg->w; mpos.h = msg->h;
		SDL_Texture *txt;
		txt = SDL_CreateTextureFromSurface(renderer, msg);	
		SDL_RenderCopy(renderer, txt, NULL, &mpos);
		SDL_FreeSurface(msg);
		SDL_DestroyTexture(txt);
		
		//drawing grids:
		int zbx = zoom_chart->getX();//zbx: Zoom chart Begins at X, other names correspondingly
		int zex = zoom_chart->getX() + zoom_chart->getSizeX();
		int zby = zoom_chart->getY();
		int zey = zoom_chart->getY() + zoom_chart->getSizeY();

		int mbx = main_chart->getX();
		int mex = main_chart->getX() + main_chart->getSizeX();
		int mby = main_chart->getY();
		int mey = main_chart->getY() + main_chart->getSizeY();
		
		for(float ygrid = -90; ygrid <= 0; ygrid += 10)
		{
			int gy = main_chart->getValueY(ygrid);
			if(gy < mby || gy > mey) continue;
			sprintf(outstr, "%d", (int)ygrid);
			msg = TTF_RenderText_Solid(font, outstr, textColor); 
			mpos.x = main_chart->getX() - msg->w; 
			mpos.y = gy - msg->h*0.5;
			mpos.w = msg->w; mpos.h = msg->h;
			SDL_Texture *txt;
			txt = SDL_CreateTextureFromSurface(renderer, msg);	
			SDL_RenderCopy(renderer, txt, NULL, &mpos);
			SDL_FreeSurface(msg);
			SDL_DestroyTexture(txt);
		}
		for(float ygrid = -90; ygrid <= 0; ygrid += 10)
		{
			int gy = zoom_chart->getValueY(ygrid);
			if(gy < zby || gy > zey) continue;
			sprintf(outstr, "%d", (int)ygrid);
			msg = TTF_RenderText_Solid(font, outstr, textColor); 
			mpos.x = zoom_chart->getX() - msg->w; 
			mpos.y = gy - msg->h*0.5;
			mpos.w = msg->w; mpos.h = msg->h;
			SDL_Texture *txt;
			txt = SDL_CreateTextureFromSurface(renderer, msg);	
			SDL_RenderCopy(renderer, txt, NULL, &mpos);
			SDL_FreeSurface(msg);
			SDL_DestroyTexture(txt);
		}
		float mbf = main_chart_freq->getV(main_chart_freq->getDataSize()-1);
		float mef = main_chart_freq->getV(1);
		for(float xgrid = 0; xgrid < 1; xgrid += 0.02)
		{
			float freq = mbf + xgrid * (mef - mbf);
			sprintf(outstr, "%.1f", freq*0.000001);
			msg = TTF_RenderText_Solid(font, outstr, textColor); 
			cpos.x = 0;
			cpos.y = 0;
			mpos.x = mbx + xgrid * (mex - mbx) - msg->h*0.5;
			mpos.y = mey + msg->w + 5;
			mpos.w = msg->w; mpos.h = msg->h;
			SDL_Texture *txt;
			txt = SDL_CreateTextureFromSurface(renderer, msg);	
			SDL_RenderCopyEx(renderer, txt, NULL, &mpos, -90, &cpos, SDL_FLIP_NONE);
			SDL_FreeSurface(msg);
			SDL_DestroyTexture(txt);
		}
		float zbf = zoom_chart_freq->getV(zoom_chart_freq->getDataSize()-1);
		float zef = zoom_chart_freq->getV(1);
		for(float xgrid = 0; xgrid < 1; xgrid += 0.05)
		{
			float freq = zbf + xgrid * (zef - zbf);
			sprintf(outstr, "%.2f", freq*0.000001);
			msg = TTF_RenderText_Solid(font, outstr, textColor); 
			cpos.x = 0;
			cpos.y = 0;
			mpos.x = zbx + xgrid * (zex - zbx) - msg->h*0.5;
			mpos.y = zey + msg->w + 5;
			mpos.w = msg->w; mpos.h = msg->h;
			SDL_Texture *txt;
			txt = SDL_CreateTextureFromSurface(renderer, msg);	
			SDL_RenderCopyEx(renderer, txt, NULL, &mpos, -90, &cpos, SDL_FLIP_NONE);
			SDL_FreeSurface(msg);
			SDL_DestroyTexture(txt);
		}
		if(debug_print) printf("text out end\n");
		
		SDL_RenderPresent(renderer);
	}
	free(drawPix);
	TTF_CloseFont( font ); 
	TTF_Quit();
	SDL_Quit();
	
	return 0;
}
