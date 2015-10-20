#ifndef ATOMS_H
#define ATOMS_H

enum tag_t {
	TAG_COOKIE,
	TAG_VALUE
};
typedef struct {
	enum tag_t tag;
	union {
		xcb_intern_atom_cookie_t cookie;
		xcb_atom_t atom;
	} u;
} intern_atom_fast_cookie_t;

xcb_atom_t intern_atom_predefined(uint16_t name_len, const char *name);
intern_atom_fast_cookie_t intern_atom_fast(xcb_connection_t *c, uint8_t only_if_exists, uint16_t name_len, const char *name);
xcb_atom_t intern_atom_fast_reply(xcb_connection_t *c, intern_atom_fast_cookie_t cookie, xcb_generic_error_t **e);

const char *get_atom_name_predefined(xcb_atom_t atom);
int get_atom_name(xcb_connection_t *c, xcb_atom_t atom, const char **namep, int *lengthp);

extern const xcb_atom_t PRIMARY;
extern const xcb_atom_t SECONDARY;
extern const xcb_atom_t ARC;
extern const xcb_atom_t ATOM;
extern const xcb_atom_t BITMAP;
extern const xcb_atom_t CARDINAL;
extern const xcb_atom_t COLORMAP;
extern const xcb_atom_t CURSOR;
extern const xcb_atom_t CUT_BUFFER0;
extern const xcb_atom_t CUT_BUFFER1;
extern const xcb_atom_t CUT_BUFFER2;
extern const xcb_atom_t CUT_BUFFER3;
extern const xcb_atom_t CUT_BUFFER4;
extern const xcb_atom_t CUT_BUFFER5;
extern const xcb_atom_t CUT_BUFFER6;
extern const xcb_atom_t CUT_BUFFER7;
extern const xcb_atom_t DRAWABLE;
extern const xcb_atom_t FONT;
extern const xcb_atom_t INTEGER;
extern const xcb_atom_t PIXMAP;
extern const xcb_atom_t POINT;
extern const xcb_atom_t RECTANGLE;
extern const xcb_atom_t RESOURCE_MANAGER;
extern const xcb_atom_t RGB_COLOR_MAP;
extern const xcb_atom_t RGB_BEST_MAP;
extern const xcb_atom_t RGB_BLUE_MAP;
extern const xcb_atom_t RGB_DEFAULT_MAP;
extern const xcb_atom_t RGB_GRAY_MAP;
extern const xcb_atom_t RGB_GREEN_MAP;
extern const xcb_atom_t RGB_RED_MAP;
extern const xcb_atom_t STRING;
extern const xcb_atom_t VISUALID;
extern const xcb_atom_t WINDOW;
extern const xcb_atom_t WM_COMMAND;
extern const xcb_atom_t WM_HINTS;
extern const xcb_atom_t WM_CLIENT_MACHINE;
extern const xcb_atom_t WM_ICON_NAME;
extern const xcb_atom_t WM_ICON_SIZE;
extern const xcb_atom_t WM_NAME;
extern const xcb_atom_t WM_NORMAL_HINTS;
extern const xcb_atom_t WM_SIZE_HINTS;
extern const xcb_atom_t WM_ZOOM_HINTS;
extern const xcb_atom_t MIN_SPACE;
extern const xcb_atom_t NORM_SPACE;
extern const xcb_atom_t MAX_SPACE;
extern const xcb_atom_t END_SPACE;
extern const xcb_atom_t SUPERSCRIPT_X;
extern const xcb_atom_t SUPERSCRIPT_Y;
extern const xcb_atom_t SUBSCRIPT_X;
extern const xcb_atom_t SUBSCRIPT_Y;
extern const xcb_atom_t UNDERLINE_POSITION;
extern const xcb_atom_t UNDERLINE_THICKNESS;
extern const xcb_atom_t STRIKEOUT_ASCENT;
extern const xcb_atom_t STRIKEOUT_DESCENT;
extern const xcb_atom_t ITALIC_ANGLE;
extern const xcb_atom_t X_HEIGHT;
extern const xcb_atom_t QUAD_WIDTH;
extern const xcb_atom_t WEIGHT;
extern const xcb_atom_t POINT_SIZE;
extern const xcb_atom_t RESOLUTION;
extern const xcb_atom_t COPYRIGHT;
extern const xcb_atom_t NOTICE;
extern const xcb_atom_t FONT_NAME;
extern const xcb_atom_t FAMILY_NAME;
extern const xcb_atom_t FULL_NAME;
extern const xcb_atom_t CAP_HEIGHT;
extern const xcb_atom_t WM_CLASS;
extern const xcb_atom_t WM_TRANSIENT_FOR;

#endif /* ATOMS_H */
