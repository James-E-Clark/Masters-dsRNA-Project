// Hive Heating Widget for Scriptable
// Displays temps in widget, triggers boost when tapped

// ========== CONFIGURATION ==========
const HA_URL = "http://100.90.18.89:8123";
const HA_TOKEN = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiIzMmFiY2NjZDAxOGE0NTFkYmIzYTY1YzIwNTJkOGQ4NCIsImlhdCI6MTc3NDk2MDAzMSwiZXhwIjoyMDkwMzIwMDMxfQ.BdHwux9VMH_wSxYdRDQx5gweb3uRAkLuSg9Ln3Vs4Yg";

const INTERNAL_TEMP_ENTITY = "sensor.hive_heating_current_temperature";
const EXTERNAL_TEMP_ENTITY = "sensor.front_temperature";
const CLIMATE_ENTITY = "climate.hive_heating_climate";
// ====================================

function roundTemp(val) {
  const num = parseFloat(val);
  return isNaN(num) ? "N/A" : num.toFixed(1);
}

async function getState(entityId) {
  const req = new Request(`${HA_URL}/api/states/${entityId}`);
  req.headers = {
    "Authorization": `Bearer ${HA_TOKEN}`,
    "Content-Type": "application/json"
  };
  try {
    const resp = await req.loadJSON();
    return resp.state;
  } catch (e) {
    return "N/A";
  }
}

async function getPresetMode() {
  const req = new Request(`${HA_URL}/api/states/${CLIMATE_ENTITY}`);
  req.headers = {
    "Authorization": `Bearer ${HA_TOKEN}`,
    "Content-Type": "application/json"
  };
  try {
    const resp = await req.loadJSON();
    return resp.attributes.preset_mode || "none";
  } catch (e) {
    return "none";
  }
}

async function triggerBoost() {
  const req = new Request(`${HA_URL}/api/services/climate/set_preset_mode`);
  req.method = "POST";
  req.headers = {
    "Authorization": `Bearer ${HA_TOKEN}`,
    "Content-Type": "application/json"
  };
  req.body = JSON.stringify({
    entity_id: CLIMATE_ENTITY,
    preset_mode: "boost"
  });
  try {
    await req.loadString();
    return true;
  } catch (e) {
    return false;
  }
}

async function cancelBoost() {
  const req = new Request(`${HA_URL}/api/services/climate/set_preset_mode`);
  req.method = "POST";
  req.headers = {
    "Authorization": `Bearer ${HA_TOKEN}`,
    "Content-Type": "application/json"
  };
  req.body = JSON.stringify({
    entity_id: CLIMATE_ENTITY,
    preset_mode: "none"
  });
  try {
    await req.loadString();
    return true;
  } catch (e) {
    return false;
  }
}

// ========== NOT IN WIDGET = TOGGLE BOOST ==========
if (!config.runsInWidget) {
  const currentMode = await getPresetMode();
  const isBoosting = currentMode === "boost";

  const success = isBoosting ? await cancelBoost() : await triggerBoost();

  const n = new Notification();
  n.title = "Heating Boost";
  if (isBoosting) {
    n.body = success ? "Boost deactivated!" : "Failed to deactivate boost";
  } else {
    n.body = success ? "Boost activated!" : "Failed to activate boost";
  }
  n.schedule();
  Script.complete();
  return;
}

// ========== IN WIDGET = SHOW TEMPS ==========
const internalTemp = roundTemp(await getState(INTERNAL_TEMP_ENTITY));
const externalTemp = roundTemp(await getState(EXTERNAL_TEMP_ENTITY));

function buildCircular() {
  const w = new ListWidget();
  w.setPadding(0, 0, 0, 0);

  const stack = w.addStack();
  stack.layoutVertically();
  stack.addSpacer(null);

  const inRow = stack.addStack();
  inRow.layoutHorizontally();
  inRow.addSpacer(null);
  const inText = inRow.addText(`🌡${internalTemp}°`);
  inText.font = Font.boldSystemFont(10);
  inText.textColor = Color.white();
  inRow.addSpacer(null);

  stack.addSpacer(1);

  const outRow = stack.addStack();
  outRow.layoutHorizontally();
  outRow.addSpacer(null);
  const outText = outRow.addText(`❄${externalTemp}°`);
  outText.font = Font.boldSystemFont(10);
  outText.textColor = Color.white();
  outRow.addSpacer(null);

  stack.addSpacer(1);

  const fireRow = stack.addStack();
  fireRow.layoutHorizontally();
  fireRow.addSpacer(null);
  const fire = fireRow.addText("🔥");
  fire.font = Font.systemFont(12);
  fireRow.addSpacer(null);

  stack.addSpacer(null);
  return w;
}

function buildRectangular() {
  const w = new ListWidget();
  w.setPadding(4, 8, 4, 8);

  const titleRow = w.addStack();
  titleRow.layoutHorizontally();
  titleRow.centerAlignContent();
  const titleIcon = titleRow.addText("🏠");
  titleIcon.font = Font.systemFont(10);
  titleRow.addSpacer(2);
  const titleText = titleRow.addText("Heating");
  titleText.font = Font.boldSystemFont(10);
  titleText.textColor = Color.white();

  w.addSpacer(2);

  const intRow = w.addStack();
  intRow.layoutHorizontally();
  intRow.centerAlignContent();
  const intIcon = intRow.addText("🌡️");
  intIcon.font = Font.systemFont(9);
  intRow.addSpacer(2);
  const intLabel = intRow.addText("In");
  intLabel.font = Font.mediumSystemFont(10);
  intLabel.textColor = Color.white();
  intRow.addSpacer(null);
  const intVal = intRow.addText(`${internalTemp}°`);
  intVal.font = Font.boldSystemFont(12);
  intVal.textColor = Color.white();

  w.addSpacer(1);

  const extRow = w.addStack();
  extRow.layoutHorizontally();
  extRow.centerAlignContent();
  const extIcon = extRow.addText("❄️");
  extIcon.font = Font.systemFont(9);
  extRow.addSpacer(2);
  const extLabel = extRow.addText("Out");
  extLabel.font = Font.mediumSystemFont(10);
  extLabel.textColor = Color.white();
  extRow.addSpacer(null);
  const extVal = extRow.addText(`${externalTemp}°`);
  extVal.font = Font.boldSystemFont(12);
  extVal.textColor = Color.white();

  return w;
}

function buildMedium() {
  const w = new ListWidget();
  w.backgroundColor = new Color("#1c1c1e");
  w.setPadding(12, 14, 12, 14);

  const title = w.addText("🏠 Heating");
  title.font = Font.boldSystemFont(13);
  title.textColor = new Color("#ff9500");

  w.addSpacer(6);

  const intStack = w.addStack();
  intStack.layoutHorizontally();
  intStack.centerAlignContent();
  const intIcon = intStack.addText("🌡️");
  intIcon.font = Font.systemFont(14);
  intStack.addSpacer(4);
  const intLabel = intStack.addText("Inside");
  intLabel.font = Font.mediumSystemFont(12);
  intLabel.textColor = new Color("#aaaaaa");
  intStack.addSpacer(null);
  const intValue = intStack.addText(`${internalTemp}°C`);
  intValue.font = Font.boldSystemFont(16);
  intValue.textColor = new Color("#ff6b6b");

  w.addSpacer(4);

  const extStack = w.addStack();
  extStack.layoutHorizontally();
  extStack.centerAlignContent();
  const extIcon = extStack.addText("❄️");
  extIcon.font = Font.systemFont(14);
  extStack.addSpacer(4);
  const extLabel = extStack.addText("Outside");
  extLabel.font = Font.mediumSystemFont(12);
  extLabel.textColor = new Color("#aaaaaa");
  extStack.addSpacer(null);
  const extValue = extStack.addText(`${externalTemp}°C`);
  extValue.font = Font.boldSystemFont(16);
  extValue.textColor = new Color("#64b5f6");

  w.addSpacer(6);

  const boostStack = w.addStack();
  boostStack.layoutHorizontally();
  boostStack.centerAlignContent();
  boostStack.backgroundColor = new Color("#ff9500", 0.2);
  boostStack.cornerRadius = 8;
  boostStack.setPadding(6, 0, 6, 0);
  boostStack.addSpacer(null);
  const boostText = boostStack.addText("🔥 TAP TO BOOST");
  boostText.font = Font.boldSystemFont(12);
  boostText.textColor = new Color("#ff9500");
  boostStack.addSpacer(null);

  return w;
}

let widget;
switch (config.widgetFamily) {
  case "accessoryCircular":
    widget = buildCircular();
    break;
  case "accessoryRectangular":
    widget = buildRectangular();
    break;
  default:
    widget = buildMedium();
    break;
}

Script.setWidget(widget);
Script.complete();
