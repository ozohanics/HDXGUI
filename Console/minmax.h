#if !defined( __MINMAX_DEFINED)
#define __MINMAX_DEFINED
#define NOMINMAX    /* for WINDEF.H */
#undef min  // make sure these aren't macros
#undef max
/*
  MFC code sometimes contains usages of min() and max() of dis-similar types
  which prevents using a normal template implementation.  We cannot
  implement min and max as macros (like Microsoft does), since parts of the
  Rogue Wave Standard Library need to #undef them.

  So we start by providing the normal template implementation and then also
  provide a global, non-template version, of min and max that take and
  return unsigned longs. The theory is that the compiler will first look at
  the non-template version and promote both params to unsigned long before
  looking at the template version and failing because of the two different
  types involved.
*/
template <class _T> inline const _T _FAR &min( const _T _FAR &__t1, const _T _FAR &__t2 )
{
	if  (__t1 < __t2)
		return __t1;
	else
		return __t2;
}

template <class _T> inline const _T _FAR &max( const _T _FAR &__t1, const _T _FAR &__t2 )
{
	if  (__t1 > __t2)
		return __t1;
	else
        return __t2;
}

inline unsigned long min (unsigned long __a, unsigned long __b)
{
	return min<unsigned long> (__a, __b);
}

inline unsigned long max (unsigned long __a, unsigned long __b)
{
	return max<unsigned long> (__a, __b);
}
#define __max       max
#define __min       min
#endif /*__MINMAX_DEFINED */

