

int invert() {
  unsigned int op_id = 0;
  unsigned int write_prop = 0;
  unsigned int index_start = 0;
  g_mu = 0.;
  // convert to even/odd geometry
  convert_lexic_to_eo(operator_list[op_id].sr0, operator_list[op_id].sr1, (spinor*) source);
  
  operator_list[op_id].inverter(op_id, index_start, write_prop);
  
  // convert back to lexicographic geometry
  convert_eo_to_lexic((spinor*) propagator, operator_list[op_id].prop0, operator_list[op_id].prop1);


  return(0);
}
