#ifndef LISTS_H
#define LISTS_H

struct list_head {
    void *next;
    void *prev;
};

/*
 * TODO: add possibility to remove item from list
 * without PITA inside the 'foreach' cycle
 */

#define init_list(x) {(x).next=&(x); (x).prev=&(x);}
#define list_empty(x) ((x).next == &(x))
/* Note that you can't so easily delete element from list while foreach'ing
   it. See del_from_list_in_foreach_and_run() macro for a typical workaround */
#define del_from_list(x) do { \
        ((struct list_head*)(x)->next)->prev=(x)->prev; \
        ((struct list_head*)(x)->prev)->next=(x)->next; \
    } while(0)

/* Safely delete and run arbitrary code then - works for foreach() only, not foreachbask()!*/
#define del_from_list_in_foreach_and_run(x, run) do { \
    struct list_head *tmp = (struct list_head*)(x)->prev; \
    ((struct list_head*)(x)->next)->prev=(x)->prev; \
    ((struct list_head*)(x)->prev)->next=(x)->next; \
    run; \
    x = (typeof(x))tmp; \
    } while(0)
/* Safely delete and run arbitrary code then - works for foreachback() only, not foreach()!*/
#define del_from_list_in_foreachback_and_run(x, run) do { \
    struct list_head *tmp = (struct list_head*)(x)->next; \
    ((struct list_head*)(x)->next)->prev=(x)->prev; \
    ((struct list_head*)(x)->prev)->next=(x)->next; \
    run; \
    x = (typeof(x))tmp; \
    } while(0)

#define add_at_pos(p,x) do { \
        (x)->next=(p)->next; \
        (x)->prev=(p); \
        (p)->next=(x); \
        (x)->next->prev=(x); \
    } while(0)

/* adds new element to the head of list, so it will be first in foreach iterator */
#define add_to_list(l,x) add_at_pos((typeof(x))&(l),(x))
/* adds it to the end of list, so it will be last in foreach cycle */
#define add_to_list_end(l,x) add_at_pos((typeof(x))list_first_item(l), (x))

/* Append the second list to the end of the first one */
#define add_list_to_list(l0, l) do {                                    \
        ((struct list_head *)(l0).prev)->next = (l).next;               \
        (l0).prev = (l).prev;                                           \
        ((struct list_head *)(l).prev)->next = ((struct list_head *)(l0).next)->prev; \
        (l).prev = &(l);                                                \
        (l).next = &(l);                                                \
    } while(0)

#define foreach(e,l) for ((e)=(typeof(e))(l).next; (e)!=(typeof(e))&(l); (e)=(typeof(e))(e)->next)
#define foreachback(e,l) for ((e)=(typeof(e))(l).prev; (e)!=(typeof(e))&(l); (e)=(typeof(e))(e)->prev)
#define free_list(l) { \
        while ((l).next != &(l)) { \
            struct list_head *a__=(struct list_head *)(l).next; \
            del_from_list(a__); free(a__);\
        } \
    }

/* Oldest item */
#define list_first_item(l) ((l).prev)
/* Newest item - last added */
#define list_last_item(l) ((l).next)

#define NULL_LIST_HEAD NULL, NULL
#define D_LIST_HEAD(x) &x, &x
#define INIT_LIST_HEAD(x) struct list_head x = { D_LIST_HEAD(x) }
#define LIST_HEAD(x) x *next; x *prev

static inline int list_length(struct list_head *list)
{
    struct list_element_str {
        LIST_HEAD(struct list_element_str);
    } *p;
    int len = 0;

    foreach(p, *list)
        len ++;

    return len;
}

#endif /* LISTS_H */
