#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "call_var_main.h"
#include "utils.h"
#include "htslib/kstring.h"
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>

// print a backtrace on a fatal signal (no debugger needed); addresses resolve with addr2line -e <binary>
static void longcalld_crash_handler(int sig) {
    void *frames[64];
    int n = backtrace(frames, 64);
    const char msg[] = "\n[longcallD] fatal signal, backtrace (resolve with addr2line -f -e bin/longcallD <addr>):\n";
    if (write(2, msg, sizeof(msg) - 1) < 0) {}
    backtrace_symbols_fd(frames, n, 2);
    signal(sig, SIG_DFL);
    raise(sig);
}

const char PROG[20] = "longcallD";
const char DESCRIP[100] = "local-haplotagging-based small and structural variant calling";
#ifndef LONGCALLD_VERSION
const char LONGCALLD_VERSION[20] = "0.0.11";
#endif
const char CONTACT[50] = "yangao07@hit.edu.cn";
int LONGCALLD_VERBOSE = 0;
char *CMD;

static int usage(void) {//main usage
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: %s: %s\n", PROG, DESCRIP);
    fprintf(stderr, "Version: %s\tContact: %s\n\n", LONGCALLD_VERSION, CONTACT); 

    fprintf(stderr, "Usage:   %s <command> [options]\n\n", PROG);

    fprintf(stderr, "Command: \n");
    fprintf(stderr, "         call          call variants from long-read BAM/CRAM, single normal sample\n");
    // fprintf(stderr, "         trio          call variants from long-read BAM/CRAM, trio normal samples\n");
    // fprintf(stderr, "         joint         joint variant calling for multiple samples\n");
    // fprintf(stderr, "         genotype      call genotype for given VCF\n");
    // fprintf(stderr, "         phase         phase given VCF based on long-read BAM\n");
    // fprintf(stderr, "         realn         refine base-level alignment\n");

    fprintf(stderr, "\n");
    return 1;
}

int main(int argc, char *argv[]) {
    signal(SIGSEGV, longcalld_crash_handler); signal(SIGBUS, longcalld_crash_handler); signal(SIGABRT, longcalld_crash_handler); signal(SIGFPE, longcalld_crash_handler);
    int i; kstring_t cmd = {0,0,0};
    kputs(argv[0], &cmd);
    for (i = 1; i < argc; ++i) {
        kputs(" ", &cmd); kputs(argv[i], &cmd);
    } CMD = strdup(cmd.s); free(cmd.s);
    
    int ret = 0;
    if (argc < 2) {
        ret = 1; usage();
    } else {
        if (strcmp(argv[1], "call") == 0) ret = call_var_main(argc-1, argv+1);
        else {
            _err_error("Unrecognized command '%s'\n", argv[1]);
            ret = 1; usage();
        }
    }
    free(CMD); return ret;
}
