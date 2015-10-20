/*
 * This work is held in copyright as an unpublished work by Silicon
 * Graphics, Inc. All rights reserved.
 */
#ifndef _ESPLOGMON_H
#define _ESPLOGMON_H

#ifdef __cplusplus
extern "C" {
#endif

int
LogEspEvent(int EvClass, int EvType, int EvVersion, char *Date, char *HostName,
        char *AppName, char *Origin, char *Msg, int RepeatCount);

#ifdef __cplusplus
}
#endif

#endif /* _ESPLOGMON_H */
