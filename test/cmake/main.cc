#ifdef CHECK_sat
bool sat_check();
#endif
#ifdef CHECK_esop
bool esop_check();
#endif

int main()
{
#ifdef CHECK_sat
  if ( !sat_check() )
    return 1;
#endif
#ifdef CHECK_esop
  if ( !esop_check() )
    return 2;
#endif
  return 0;
}
