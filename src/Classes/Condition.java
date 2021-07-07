package Classes;

public class Condition extends Item{
    //public Condition() {
  //  }
// CONSTRUCTOR LLENO
  //  public Condition(float z, int node4) {
    //    super(z, node4);
    //}

    //condition requiere un constructor lleno no vacio

    //Metodo que nos ayudara a crear las listas de condiciones
    public static Condition[] createConditions(int n){
        Condition[] list = new Condition[n];
        for (int i = 0; i < n; i++) {
            list[i] = new Condition();
        }
        return list;
    }

    //la cantidad de parametros que se le manda a set value son mas valores
    @Override
    public void setValues(int a,float b,float c,float d,int e,int f,int g, int h, float i) {
        node1 = e;
        value = i;
    }



}
