/*
 * Copyright (c) 2010 Ken McDonell.  All Rights Reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
 */

#ifndef _IMPORT_H
#define _IMPORT_H

#ifdef __cplusplus
extern "C" {
#endif

// core libpcp_import API routines
//
extern int pmiStart(const char *, int);
extern int pmiUseContext(int);
extern int pmiEnd(void);
extern int pmiSetHostname(const char *);
extern int pmiSetTimezone(const char *);
extern int pmiAddMetric(const char *, pmID, int, pmInDom, int, pmUnits);
extern int pmiAddInstance(pmInDom, const char *, int);
extern int pmiPutValue(const char *, const char *, const char *);
extern int pmiGetHandle(const char *, const char *);
extern int pmiPutValueHandle(int, const char *);
extern int pmiWrite(int, int);
extern int pmiPutResult(const pmResult *);

// helper functions
//
extern void pmiDump(void);
extern pmUnits pmiUnits(int, int, int, int, int, int);
extern const char *pmiErrStr(int);

// libpcp_import error codes
//

#define PMI_ERR_BASE 20000
#define PMI_ERR_DUPMETRICNAME	(-PMI_ERR_BASE-1) /* Metric name already defined */
#define PMI_ERR_DUPMETRICID	(-PMI_ERR_BASE-2) /* Metric pmID already defined */
#define PMI_ERR_DUPINSTNAME	(-PMI_ERR_BASE-3) /* External instance name already defined */
#define PMI_ERR_DUPINSTID	(-PMI_ERR_BASE-4) /* Internal instance identifer already defined */
#define PMI_ERR_INSTNOTNULL	(-PMI_ERR_BASE-5) /* Non-null instance expected for a singular metric */
#define PMI_ERR_INSTNULL	(-PMI_ERR_BASE-6) /* Null instance not allowed for a non-singular metric */
#define PMI_ERR_BADHANDLE	(-PMI_ERR_BASE-7) /* Illegal handle */
#define PMI_ERR_DUPVALUE	(-PMI_ERR_BASE-8) /* Value already assigned for singular metric */
#define PMI_ERR_BADTYPE		(-PMI_ERR_BASE-9) /* Illegal metric type */
#define PMI_ERR_BADSEM		(-PMI_ERR_BASE-10) /* Illegal metric semantics */
#define PMI_ERR_NODATA		(-PMI_ERR_BASE-11) /* No data to output */

#ifdef __cplusplus
}
#endif

#endif /* _IMPORT_H */
