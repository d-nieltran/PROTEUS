import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import ProteusVS from "./proteus-v4.jsx";

createRoot(document.getElementById("root")).render(
  <StrictMode>
    <ProteusVS />
  </StrictMode>
);
