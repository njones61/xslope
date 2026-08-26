// Collapsed sidebar for the ReadTheDocs theme: each top-level category (rendered
// as a caption followed by a list) shows its pages only while one of them is the
// current page. Clicking a caption toggles it.
document.addEventListener("DOMContentLoaded", function () {
  var captions = document.querySelectorAll(".wy-menu-vertical p.caption");
  captions.forEach(function (cap) {
    var list = cap.nextElementSibling;
    if (!list || list.tagName !== "UL") return;
    var open = list.querySelector("li.current") !== null;
    cap.classList.add("nav-collapsible");
    if (!open) { list.style.display = "none"; cap.classList.add("nav-closed"); }
    cap.addEventListener("click", function () {
      var hidden = list.style.display === "none";
      list.style.display = hidden ? "" : "none";
      cap.classList.toggle("nav-closed", !hidden);
    });
  });
});
