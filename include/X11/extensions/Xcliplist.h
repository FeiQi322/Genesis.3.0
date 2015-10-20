#ifndef _XCLIPLIST_H_
#define _XCLIPLIST_H_

#include <X11/Xfuncproto.h>

#define X_XClipListQueryVersion		0
#define X_XGetClipList			1

#define XClipListNumberEvents		0

#define XClipListClientNotLocal		0
#define XClipListOperationNotSupported	1
#define XClipListNumberErrors		(XClipListOperationNotSupported + 1)

_XFUNCPROTOBEGIN

#ifndef _XCLIPLIST_SERVER_

typedef struct _BoxRec {
    short x1, y1, x2, y2;
} BoxRec, *BoxPtr;


Bool XClipListQueryExtension(
    Display*		/* dpy */,
    int*		/* event_base */,
    int*		/* error_base */
);

Bool XClipListQueryVersion(
    Display*		/* dpy */,
    int*		/* majorVersion */,
    int*		/* minorVersion */,
    int*		/* patchVersion */
);

BoxPtr XGetClipList(
    Display*		/* dpy */,
    XID			/* Window ID */,
    int*		/* number of clip regions */
   );

#endif

_XFUNCPROTOEND

#endif /* _XCLIPLIST_H_ */

