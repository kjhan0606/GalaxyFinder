#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <unistd.h>
#include <sys/time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct StageTimer {
	double wall0;
	double cpu0;
} StageTimer;

static inline double stage_wall(void) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec * 1e-6;
}

/* Sum of user+system time over every thread of this process. */
static inline double stage_cpu(void) {
	DIR *dir = opendir("/proc/self/task");
	struct dirent *ent;
	long ticks;
	double sec = 0;
	if(!dir) return 0;
	ticks = sysconf(_SC_CLK_TCK);
	if(ticks <= 0) ticks = 100;
	while((ent = readdir(dir)) != NULL){
		char path[96], buf[1024], *paren;
		FILE *fp;
		unsigned long ut = 0, st = 0;
		if(ent->d_name[0] == '.') continue;
		snprintf(path, sizeof path, "/proc/self/task/%s/stat", ent->d_name);
		fp = fopen(path, "r");
		if(!fp) continue;
		if(!fgets(buf, sizeof buf, fp)){
			fclose(fp);
			continue;
		}
		fclose(fp);
		paren = strrchr(buf, ')');
		if(!paren) continue;
		if(sscanf(paren + 2,
				"%*c %*d %*d %*d %*d %*d %*u %*u %*u %*u %*u %lu %lu",
				&ut, &st) >= 2)
			sec += (double)(ut + st) / (double)ticks;
	}
	closedir(dir);
	return sec;
}

static inline void stage_begin(StageTimer *t) {
	t->wall0 = stage_wall();
	t->cpu0 = stage_cpu();
}

static inline void stage_end(const StageTimer *t, const char *name) {
	double wall = stage_wall() - t->wall0;
	double cpu = stage_cpu() - t->cpu0;
	int threads = 1;
#ifdef _OPENMP
	threads = omp_get_max_threads();
#endif
	{
		double eff = (wall > 0.0 && threads > 0) ? cpu / (wall * (double)threads) : 0.0;
		fprintf(stderr, "TIMER %s wall=%.3f cpu=%.3f threads=%d eff=%.3f\n",
				name, wall, cpu, threads, eff);
		fflush(stderr);
	}
}
