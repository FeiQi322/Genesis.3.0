/*
 * built from impl.h on Tue Mar 27 01:22:22 UTC 2012 by mk_pmdbg
 */

typedef struct {
    char	*name;
    int		bit;
} debug_map_t;

static debug_map_t debug_map[] = {
    { "AF",	DBG_TRACE_AF },
    { "APPL0",	DBG_TRACE_APPL0 },
    { "APPL1",	DBG_TRACE_APPL1 },
    { "APPL2",	DBG_TRACE_APPL2 },
    { "CONFIG",	DBG_TRACE_CONFIG },
    { "CONTEXT",	DBG_TRACE_CONTEXT },
    { "DERIVE",	DBG_TRACE_DERIVE },
    { "FETCH",	DBG_TRACE_FETCH },
    { "INDOM",	DBG_TRACE_INDOM },
    { "INTERP",	DBG_TRACE_INTERP },
    { "LIBPMDA",	DBG_TRACE_LIBPMDA },
    { "LOG",	DBG_TRACE_LOG },
    { "LOGMETA",	DBG_TRACE_LOGMETA },
    { "LOOP",	DBG_TRACE_LOOP },
    { "OPTFETCH",	DBG_TRACE_OPTFETCH },
    { "PDU",	DBG_TRACE_PDU },
    { "PDUBUF",	DBG_TRACE_PDUBUF },
    { "PMC",	DBG_TRACE_PMC },
    { "PMNS",	DBG_TRACE_PMNS },
    { "PROFILE",	DBG_TRACE_PROFILE },
    { "TIMECONTROL",	DBG_TRACE_TIMECONTROL },
    { "VALUE",	DBG_TRACE_VALUE },
};

static int num_debug = sizeof(debug_map) / sizeof(debug_map[0]);
