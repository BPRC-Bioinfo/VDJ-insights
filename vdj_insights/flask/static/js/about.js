document.addEventListener("DOMContentLoaded", function () {
    const copyButtons = document.querySelectorAll(".copy-button");

    copyButtons.forEach(function (button) {
        button.addEventListener("click", function () {
            const targetId = this.getAttribute("data-target");
            const codeElement = document.getElementById(targetId);
            const code = codeElement.textContent;

            // Copy the code to clipboard
            navigator.clipboard.writeText(code).then(() => {
                // Success feedback
                const iconElement = this.querySelector("i");
                if (iconElement) {
                    iconElement.className = "fa-solid fa-check icon-active";
                }

                // Show the tooltip
                this.classList.add('show-tooltip');

                // Reset after 2 seconds
                setTimeout(() => {
                    if (iconElement) {
                        iconElement.className = "fa-regular fa-clone";
                        iconElement.classList.remove("icon-active");
                    }
                    this.classList.remove('show-tooltip');
                }, 2000);
            }, () => {
                // Error handling if needed
                console.error("Failed to copy text.");
            });
        });
    });
});
