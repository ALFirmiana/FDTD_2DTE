#include "field.h"

class simu
{
  public:
    void evol(Field *field);

  private:
    void evol_inner(int i, int j, Field *field);
    void evol_boundary();
    void evol_inject();
    Field *field1, field2;
};
