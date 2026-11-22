// "NEW" pills in the sidebar for what the latest release added. Update the two
// lists at each release. Each entry is the end of a page's URL path; a sidebar link
// is matched on the address it resolves to, so relative links and the current page
// (which the theme links as "#") match alike.
var NEW_PAGES = [                 // a pill after the page's own sidebar entry
  "/tutorials/lem14_pseudostatic/",
  "/tutorials/lem15_yield_newmark/",
];
var NEW_SECTIONS = [              // a pill after the heading of the section that holds the page
  "/seismic/overview/",
];

document.addEventListener("DOMContentLoaded", function () {
  function path(a) {
    var href = a.getAttribute("href");
    if (href === null) return "";
    if (href === "#") return window.location.pathname;
    return new URL(href, window.location.href).pathname;
  }
  function listed(a, list) {
    var p = path(a);
    return list.some(function (end) { return p.endsWith(end); });
  }
  document.querySelectorAll(".wy-menu-vertical a").forEach(function (a) {
    if (listed(a, NEW_PAGES)) a.classList.add("nav-new");
  });
  document.querySelectorAll(".wy-menu-vertical p.caption").forEach(function (cap) {
    var list = cap.nextElementSibling;
    if (!list || list.tagName !== "UL") return;
    var hit = Array.prototype.some.call(list.querySelectorAll("a"), function (a) {
      return listed(a, NEW_SECTIONS);
    });
    var text = cap.querySelector(".caption-text");
    if (hit && text) text.classList.add("nav-new");
  });
});
