/* Per-particle marks for the halo being identified.
 *
 * One byte, set and cleared by the functions below. The bits are the
 * historical values: code that dumped this byte still reads the same mask.
 * Thread-private visit marks live in a separate byte per OpenMP thread.
 */
#pragma once

#define PARTICLE_IS_PEAK      0x01
#define PARTICLE_IS_VISITED   0x02
#define PARTICLE_IS_CORE      0x04
#define PARTICLE_IS_SHELL     0x08
#define PARTICLE_IS_BOUND     0x10
#define PARTICLE_IS_REMAINING 0x20
#define THREAD_IS_VISITED     0x01

#define CORE_MARK_ENCLOSED  0x01
#define CORE_MARK_CONFIRMED 0x02
#define CORE_MARK_REAL      0x04

/* Zero, used by the old "IS_*(i) == NOT" tests. */
#define NOT 0

typedef struct ParticleState {
    float density;
    int galaxy_id;
    unsigned char marks;
    unsigned char thread_marks[MAXTHREADS];
} ParticleState;

/* Core marks live on Coretype.marks. */
static inline void reset_core_marks(Coretype *cores, int i) { cores[i].marks = 0; }
static inline void mark_core_enclosed(Coretype *cores, int i)  { cores[i].marks |= CORE_MARK_ENCLOSED; }
static inline void clear_core_enclosed(Coretype *cores, int i) { cores[i].marks &= (unsigned char)~CORE_MARK_ENCLOSED; }
static inline int  core_is_enclosed(Coretype *cores, int i)    { return cores[i].marks & CORE_MARK_ENCLOSED; }
static inline void toggle_core_enclosed(Coretype *cores, int i){ cores[i].marks ^= CORE_MARK_ENCLOSED; }
static inline void mark_real_core(Coretype *cores, int i)      { cores[i].marks |= CORE_MARK_REAL; }
static inline void clear_real_core(Coretype *cores, int i)     { cores[i].marks &= (unsigned char)~CORE_MARK_REAL; }
static inline int  core_is_real(Coretype *cores, int i)        { return cores[i].marks & CORE_MARK_REAL; }
static inline void toggle_real_core(Coretype *cores, int i)    { cores[i].marks ^= CORE_MARK_REAL; }

static inline void reset_order_marks(Coresortdentype *row)      { row->flag = 0; }
static inline void mark_order_enclosed(Coresortdentype *row)    { row->flag |= CORE_MARK_ENCLOSED; }
static inline void clear_order_enclosed(Coresortdentype *row)   { row->flag &= (unsigned char)~CORE_MARK_ENCLOSED; }
static inline int  order_is_enclosed(Coresortdentype *row)      { return row->flag & CORE_MARK_ENCLOSED; }
static inline void toggle_order_enclosed(Coresortdentype *row)  { row->flag ^= CORE_MARK_ENCLOSED; }
static inline void mark_order_confirmed(Coresortdentype *row)   { row->flag |= CORE_MARK_CONFIRMED; }
static inline void clear_order_confirmed(Coresortdentype *row)  { row->flag &= (unsigned char)~CORE_MARK_CONFIRMED; }
static inline int  order_is_confirmed(Coresortdentype *row)     { return row->flag & CORE_MARK_CONFIRMED; }
static inline void toggle_order_confirmed(Coresortdentype *row) { row->flag ^= CORE_MARK_CONFIRMED; }
