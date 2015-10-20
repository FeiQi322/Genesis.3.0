#ifndef _XCLIPLISTSTR_H_
#define _XCLIPLISTSTR_H_

#include "Xcliplist.h"

#define XCLIPLISTNAME "XClipList"

/*
 * Version History:
 *    1.0.0: Original
 */
#define XCLIPLIST_MAJOR_VERSION	1
#define XCLIPLIST_MINOR_VERSION	0
#define XCLIPLIST_PATCH_VERSION	0

typedef struct _XClipListQueryVersion {
    CARD8	reqType;		/* always DHReqCode */
    CARD8	dhReqType;		/* always X_DHQueryVersion */
    CARD16	length B16;
} xXClipListQueryVersionReq;
#define sz_xXClipListQueryVersionReq	4

typedef struct {
    BYTE	type;			/* X_Reply */
    BOOL	pad1;
    CARD16	sequenceNumber B16;
    CARD32	length B32;
    CARD16	majorVersion B16;	/* major version of DH protocol */
    CARD16	minorVersion B16;	/* minor version of DH protocol */
    CARD32	patchVersion B32;       /* patch version of DH protocol */
    CARD32	pad3 B32;
    CARD32	pad4 B32;
    CARD32	pad5 B32;
    CARD32	pad6 B32;
} xXClipListQueryVersionReply;
#define sz_xXClipListQueryVersionReply	32


typedef struct _XGetClipList {
    CARD8	reqType;		/* always DHReqCode */
    CARD8	dhReqType;		/* always X_DHSwapBuffers */
    CARD16	length B16;
    CARD32	windowid B32;
} xXGetClipListReq;
#define sz_xXGetClipListReq	8

typedef struct {
    BYTE	type;			/* X_Reply */
    BOOL	pad1;
    CARD16	sequenceNumber B16;
    CARD32	length B32;
    CARD32	num_entries B32;
    CARD32	pad2 B32;
    CARD32	pad3 B32;
    CARD32	pad4 B32;
    CARD32	pad5 B32;
    CARD32	pad6 B32;
} xXGetClipListReply;

#define sz_xXGetClipListReply	32

#endif /* _XCLIPLISTSTR_H_ */
