#include <math.h>

#include <sys/ipc.h>
#include <sys/shm.h>

#include <stdio.h>
#include <stdint.h>

#include "csvReader.h"

typedef struct sRecord
{
	float freq;
	float pwr;
	float band1;
	float band2;
}sRecord;

sRecord *records = NULL;
int records_count = 0;

void load_file(char* fname)
{
//	printf("loading: %s\n", fname);
	csvReader *R = new csvReader(fname);
	if(records_count > 0)
	{
		delete records;
		records_count = 0;
	}
	records_count = R->getLinesCount();
	records = new sRecord[records_count];
	for(int l = 0; l < records_count; l++)
	{
		records[l].freq = R->readDouble();
		records[l].pwr = R->readDouble();
		records[l].band1 = R->readDouble();
		records[l].band2 = R->readDouble();
		if(l < 10)
		{
			printf("%g %g %g %g\n", records[l].freq, records[l].pwr, records[l].band1, records[l].band2);
		}
	}
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
    }

    if ((shared_memory = (uint8_t*)shmat(shmid, NULL, 0)) == (uint8_t *) -1) {
        printf("shmat error!\n");
    }
}

int cur_file_id = 590;

void fill_scan_data()
{
	float *f_shm = (float*)shared_memory;
	int *i_shm = (int*)shared_memory;
//	float start_freq = f_shm[1];
//	float end_freq = f_shm[2];
//	float gain_mod = f_shm[3];
	i_shm[4] = records_count;
	
	for(int r = 0; r < records_count; r++)
	{
		f_shm[5 + r*2] = records[r].freq;
		f_shm[6 + r*2] = records[r].band1;
	}

	i_shm[0] = cur_file_id;
}

void shared_mem_controller()
{
	if(!shared_memory) shared_mem_init();
}

int main()
{
	shared_mem_init();
	int minute = 0;
	int second = 0;
	while(1)
	{
		char fnm[512];
		records_count = 0;
		int skip_count = 0;
		while(records_count < 1)
		{
			sprintf(fnm, "../new_logs/full_range_logs/signal_00_%02d_%02d_%d.000000_%d.000000.txt", minute, second, cur_file_id, cur_file_id+20);
			load_file(fnm);
			if(records_count < 1)
			{
				second++;
				skip_count++;
				if(second > 59)
				{
					second = 0;
					minute++;
				}
			}
			if(skip_count > 100) break;
		}
		if(records_count < 1) { minute = 1111; break; };
		printf("file %s, records: %d\n", fnm, records_count);
		if(records_count > 0)
		{
			fill_scan_data();
		}
		cur_file_id += 4;
		usleep(10000);
	}
	return 0;
}