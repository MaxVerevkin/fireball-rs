function to_degrees(a) {
  return a / Math.PI * 180.0;
}
function vec3_sub(a, b) {
  return {
    x: a.x - b.x,
    y: a.y - b.y,
    z: a.z - b.z
  };
}
function vec3_scal_mul(a, num) {
  return {
    x: a.x * num,
    y: a.y * num,
    z: a.z * num
  };
}
function vec3_norm_sq(a) {
  return a.x * a.x + a.y * a.y + a.z * a.z;
}
function vec3_norm(a) {
  return Math.sqrt(vec3_norm_sq(a));
}
function vec3_normalize(a) {
  return vec3_scal_mul(a, 1.0 / vec3_norm(a));
}
function vec3_dot(a, b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
function vec3_cross(a, b) {
  return {
    x: a.y * b.z - a.z * b.y,
    y: a.z * b.x - a.x * b.z,
    z: a.x * b.y - a.y * b.x
  };
}
function vec3_angle_degrees(a, b) {
  return to_degrees(Math.acos(Math.max(Math.min(vec3_dot(a, b), 1.0), -1.0)));
}
function open_db() {
  const input = document.createElement("input");
  input.type = "file";
  input.accept = ".json";
  input.onchange = (e)=>{
    const reader = new FileReader();
    // @ts-ignore
    reader.readAsText(e.target.files[0], "UTF-8");
    // @ts-ignore
    reader.onload = (readerEvent)=>process_data(JSON.parse(readerEvent.target.result));
  };
  input.click();
}
function extract_all_params(data) {
  const all_params = [];
  for (const event of Object.values(data.events))for (const param of Object.keys(event.runs))if (!all_params.includes(param)) all_params.push(param);
  return all_params;
}
function extract_all_stages(data, params) {
  let stage_names = [];
  for (const [_, event] of Object.entries(data.events)){
    let counts = {};
    let run_stages = [];
    for (const stage of event.runs[params]){
      const count = counts[stage.title] ?? 0;
      run_stages.push([
        stage.title,
        count
      ]);
      counts[stage.title] = count + 1;
    }
    for (const stage of run_stages){
      if (stage_names.find((a)=>a[0] == stage[0] && a[1] == stage[1]) === undefined) {
        stage_names.push(stage);
      }
    }
  }
  return stage_names;
}
function sigma_to_class_name(error, sigma) {
  if (error / sigma > 3) return "more-than-three-sigmas";
  if (error / sigma > 2) return "two-to-three-sigmas";
  return "less-than-three-sigmas";
}
function process_data(data, params_) {
  const all_params = extract_all_params(data);
  if (all_params.length == 0) {
    window.alert("could not find any params");
    return;
  }
  const params = params_ ?? all_params[0];
  const stage_names = extract_all_stages(data, params);
  const events = Object.entries(data.events).filter((a)=>params in a[1].runs).sort((a, b)=>b[1].observers_count - a[1].observers_count);
  const table = document.createElement("table");
  const table_header = document.createElement("tr");
  const table_header2 = document.createElement("tr");
  table_header.innerHTML = `<th rowspan="2">Event (${events.length})</th> <th rowspan="2">#</th> `;
  table.appendChild(table_header);
  table.appendChild(table_header2);
  const COL_SPAN = 4;
  for (const [stage_name, stage_i] of stage_names){
    let title = stage_name;
    if (stage_i != 0) title += ` (${stage_i + 1})`;
    const cell = document.createElement("th");
    cell.innerText = title;
    cell.classList.add("bold-left-border");
    cell.colSpan = COL_SPAN;
    table_header.appendChild(cell);
    [
      "Distance, km",
      "Angle, °",
      "El. Angle, °",
      "Info"
    ].forEach((name, i)=>{
      const cell = document.createElement("th");
      cell.classList.add("nowrap");
      if (i == 0) cell.classList.add("bold-left-border");
      cell.innerHTML = name;
      table_header2.appendChild(cell);
    });
  }
  for (const [event_name, event] of events){
    const title_cell = document.createElement("td");
    const number_cell = document.createElement("td");
    title_cell.innerText = event_name;
    title_cell.classList.add("nowrap");
    number_cell.innerText = event.observers_count.toString();
    number_cell.style.cssText = "text-align: right";
    const row = document.createElement("tr");
    row.appendChild(title_cell);
    row.appendChild(number_cell);
    for (const [stage_name, stage_i] of stage_names){
      const run = event.runs[params].filter((x)=>x.title == stage_name)[stage_i];
      if (run == undefined) {
        const na_cell = document.createElement("td");
        na_cell.colSpan = COL_SPAN;
        na_cell.innerText = "N/A";
        na_cell.style.cssText = "text-align: center";
        na_cell.classList.add("bold-left-border");
        row.appendChild(na_cell);
        continue;
      }
      const distance_cell = document.createElement("td");
      const distance_value = vec3_norm(vec3_sub(event.answer.point, run.traj.point)) * 1e-3;
      distance_cell.innerHTML = distance_value.toFixed(0);
      if (run.sigmas != null) {
        const sigma_value = vec3_norm(run.sigmas) * 1e-3;
        distance_cell.innerHTML += `&pm;${sigma_value.toFixed(0)}`;
        distance_cell.classList.add(sigma_to_class_name(distance_value, sigma_value));
      }
      distance_cell.style.cssText = "text-align: right";
      distance_cell.classList.add("bold-left-border");
      row.appendChild(distance_cell);
      if (event.answer.dir == null) {
        const na_cell = document.createElement("td");
        na_cell.colSpan = 2;
        na_cell.innerText = "N/A";
        na_cell.style.cssText = "text-align: center";
        row.appendChild(na_cell);
      } else {
        const angle_value = vec3_angle_degrees(event.answer.dir, run.traj.direction);
        const el_angle_value = Math.abs(to_degrees(Math.asin(vec3_dot(vec3_cross(event.answer.dir, vec3_normalize(vec3_sub(event.avg_observer, event.answer.point))), run.traj.direction))));
        const angle_cell = document.createElement("td");
        const el_angle_cell = document.createElement("td");
        angle_cell.innerHTML = angle_value.toFixed(0);
        if (run.sigmas != null) {
          const sigma_value = to_degrees(run.sigmas.v_angle);
          angle_cell.innerHTML += `&pm;${sigma_value.toFixed(0)}`;
          angle_cell.classList.add(sigma_to_class_name(angle_value, sigma_value));
          el_angle_cell.classList.add(sigma_to_class_name(el_angle_value, sigma_value));
        }
        angle_cell.style.cssText = "text-align: right";
        row.appendChild(angle_cell);
        el_angle_cell.innerHTML = el_angle_value.toFixed(0);
        el_angle_cell.style.cssText = "text-align: right";
        row.appendChild(el_angle_cell);
      }
      if (run.extra_info == null) {
        const cell = document.createElement("td");
        cell.innerText = "N/A";
        row.appendChild(cell);
      } else {
        const cell = document.createElement("td");
        let m = run.extra_info.match(/flipped (\d+) descent angles/);
        if (m != null) {
          let n = Number(m[1]);
          let p = n / event.observers_count * 100;
          cell.innerText = `${n} (${p.toFixed(0)}%)`;
          cell.style.cssText = "text-align: right";
        } else {
          cell.innerText = run.extra_info;
        }
        row.appendChild(cell);
      }
    }
    table.appendChild(row);
  }
  const select_label = document.createElement("label");
  select_label.htmlFor = "params_select";
  select_label.innerText = "Params";
  const select = document.createElement("select");
  select.id = "params_select";
  for (const params of all_params){
    const element = document.createElement("option");
    element.value = params;
    element.innerText = params;
    select.appendChild(element);
  }
  // @ts-ignore
  select.onchange = (e)=>process_data(data, e.target.value);
  select.value = params;
  const output = document.getElementById("output");
  output.innerText = "";
  output.appendChild(select_label);
  output.appendChild(select);
  output.appendChild(document.createElement("hr"));
  output.appendChild(table);
}
