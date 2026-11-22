// "NEW" pills in the sidebar for what the latest release added. Update the two
// lists at each release. Each entry is the end of a page's URL path; a sidebar link
// is matched on the address it resolves to, so relative links and the current page
// (which the theme links as "#") match alike.
var NEW_PAGES = [       // a pill on the page's own entry and on every heading above it
  "/tutorials/lem14_pseudostatic/",
  "/tutorials/lem15_yield_newmark/",
];
var NEW_SECTIONS = [    // a pill on the headings above the page, not on the page itself
  "/seismic/overview/",
];

document.addEventListener("DOMContentLoaded", function () {
  var menu = document.querySelector(".wy-menu-vertical");
  if (!menu) return;
  function path(a) {
    var href = a.getAttribute("href");
    if (href === null) return "";
    if (href === "#") return window.location.pathname;
    return new URL(href, window.location.href).pathname;
  }
  function listed(a, list) {
    var p = path(a);
    return p !== "" && list.some(function (end) { return p.endsWith(end); });
  }
  function holds(list, ul) {
    return Array.prototype.some.call(ul.querySelectorAll("a"), function (a) {
      return listed(a, list);
    });
  }
  var all = NEW_PAGES.concat(NEW_SECTIONS);
  menu.querySelectorAll("a").forEach(function (a) {
    if (listed(a, NEW_PAGES)) a.classList.add("nav-new");
    // A heading inside a section (an entry with no page of its own) holds a list.
    var sub = a.nextElementSibling;
    if (a.getAttribute("href") === null && sub && sub.tagName === "UL" && holds(all, sub)) {
      a.classList.add("nav-new");
    }
  });
  menu.querySelectorAll("p.caption").forEach(function (cap) {
    var list = cap.nextElementSibling;
    var text = cap.querySelector(".caption-text");
    if (list && list.tagName === "UL" && text && holds(all, list)) text.classList.add("nav-new");
  });
});
